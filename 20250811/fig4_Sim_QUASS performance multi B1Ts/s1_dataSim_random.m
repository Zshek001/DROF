% 202500808
% compare CEST fitting data with different Ts/B1
%   (1) BME simulation at given Ts/B1
%   (2) CEST fitting w/wo QUASS
%   (3) linear regression
%   (4) save as .mat (one B1 for one file)
%       store zspec, sensitivity & R2, contribution

clear

% ----- Cross-platform path handling -----
here = fileparts(mfilename('fullpath'));  % folder of this script
if isempty(here)                          % e.g., if run as a Live Script
    here = pwd;
end

% Add toolbox folder (.. / toolbox)
% addpath ..\toolbox\
% addpath ..\toolbox\fitting_main_par\
toolboxDir = fullfile(here, '..', 'toolbox');
fitting_main_parDir = fullfile(here, '..', 'toolbox', 'fitting_main_par');
if exist(toolboxDir, 'dir')
    addpath(toolboxDir);
    addpath(fitting_main_parDir);
else
    warning('Toolbox folder not found: %s', toolboxDir);
end


Npixel = 5; % Zspec number
multistart_N = 1; % times for MultiStart
rng default % for reproducibility

QUASSpro_methods = @QUASSprocess_par;
fitting_methods = {@fitting_MTR_par,@fitting_LDA_par, @fitting_DMPLF_par};

% elapsedtime: ~39min for N=1000, nB1=13, nTs=12
B1List = linspace(0.3, 1.5, 5);
TsList = [10,linspace(0.8, 3.0, 5)];
outputfolder = fullfile(here, '.', 'nB1-nTs');

% delete(gcp('nocreate'));
% parpool(16);

%% check ouput folder
cnt = 1;
dirTemp = outputfolder;
while(1)
    if exist(dirTemp,'dir')
        dirTemp = outputfolder + "_" + num2str(cnt);
        cnt = cnt + 1;
    else
        break;
    end
end
outputfolder = dirTemp;

if ~exist("outputfolder",'dir')
    mkdir(outputfolder)
end

%% load freq offs and initialization
fileID = fopen("user17_77p.txt", 'r'); % included M0
data = cell2mat(textscan(fileID, '%f'));
fclose(fileID);
offs = reshape(data,[],1);
offs_trunc = offs(5:end);
nf_trunc = length(offs_trunc);

nB1 = length(B1List);
nTs = length(TsList);
TdList = TsList;
B1_Ts_ListPath = fullfile(outputfolder, 'B1_Ts_List.mat');
save(B1_Ts_ListPath,'B1List','TsList','TdList');

for idxB1 = 1:nB1
    fprintf("processing B1="+num2str(B1List(idxB1),3)+"uT: " );

    zspec_vec_givenB1 = zeros(nf_trunc,Npixel,nTs);
    zspecQUASS_vec_givenB1 = zeros(nf_trunc,Npixel,nTs);
    contri_amide_all_givenB1 = zeros(4,6,nTs); % contri from 4pools, 6 algorithm (MTR,LDA,DMPLF,Q-MTR,Q-LDA,Q-DMPLF)
    contri_guan_all_givenB1 = zeros(4,6,nTs);
    contri_noe_all_givenB1 = zeros(4,6,nTs);
    sensR2_amide_all_givenB1 = zeros(2,6,nTs); % slope and R2
    sensR2_guan_all_givenB1 = zeros(2,6,nTs);
    sensR2_noe_all_givenB1 = zeros(2,6,nTs);
    diffss_meanList_all_givenB1 = zeros(3,6,nTs); % 4pools, 6 algo

    for idxTs = 1:nTs
        backNum = fprintf(num2str(idxTs)+"/"+num2str(nTs));
        SATpara = [B1List(idxB1), TsList(idxTs), TsList(idxTs)]; % [B1(uT), tsat(s), tdelay(s)]
        
        %% simulation with M0
        [zspecList,paraList] = func1_dataSim_givenSATpara_par(SATpara, offs, Npixel); 
            %   zspecList: [nf_trunc,Npixel], truncated M0 out
            %   paraList: [4,Npixel], concentration of 4 pools
        
        %% QUASS
        SATpara = reshape(SATpara,[],1);
        
        MTf_mM = paraList(1,:);
        amidef_mM = paraList(2,:); % mM
        guanf_mM = paraList(3,:);
        noef_mM = paraList(4,:);
        
        % construct predictor variable matrix for linear regression, [nPixel,nvar]
        predictVar = [guanf_mM(:),amidef_mM(:),MTf_mM(:),noef_mM(:)];
        
        % zspec w/wo QUASS
        [nf,Npixel] = size(zspecList);
        zspec_vec = zspecList; % [nf, nzspec]
        
        % gamma_hz = 267.5153/2/pi; % for protons, in 1e6 Hz/T
        % offs_hz = offs_trunc*3*gamma_hz;
        R1w = 1; % Hz, fixed in BME simulation setting
        [zspecQUASS_vec, R1rhoList] = QUASSpro_methods(zspec_vec, offs_trunc, SATpara, R1w, 0);
            
        %% CEST fitting/analysis
        contri_amide_all = zeros(4,6); % contribution from 4pools, 6 algorithm (MTR,LDA,DMPLF,Q-MTR,Q-LDA,Q-DMPLF)
        contri_guan_all = zeros(4,6);
        contri_noe_all = zeros(4,6);
        sensR2_amide_all = zeros(2,6); % slope and R2
        sensR2_guan_all = zeros(2,6);
        sensR2_noe_all = zeros(2,6);
        if idxTs == 1
            % 1st Ts for steady-state (Ts->infty)
            ssContrast_all = zeros(3,6,Npixel); % 3pools, 6 algo
        end
        
        for idxalgo = 1:length(fitting_methods)
            method = fitting_methods{idxalgo};
            for idxprepro = 1:2 % QUASS
                if idxprepro == 1
                    zspecTemp_vec = zspec_vec;
                else
                    zspecTemp_vec = zspecQUASS_vec;
                end
                
                [Zamide_DZ_vec, Zguan_DZ_vec, Znoe_DZ_vec] = method(offs_trunc, zspecTemp_vec, multistart_N, 0);
            
                % construct response variable vector, [nPixel,1]
                responVar_guan = Zguan_DZ_vec(:);
                responVar_amide = Zamide_DZ_vec(:);
                responVar_noe = Znoe_DZ_vec(:);
            
                % linear regression
                model_4var_guan = fitlm(predictVar, responVar_guan,'ResponseVar','metric_guan','PredictorVars',{'guanf','amidef','MTf','NOEf'});
                model_4var_amide = fitlm(predictVar, responVar_amide,'ResponseVar','metric_amide','PredictorVars',{'guanf','amidef','MTf','NOEf'});
                model_4var_noe = fitlm(predictVar, responVar_noe,'ResponseVar','metric_noe','PredictorVars',{'guanf','amidef','MTf','NOEf'});
            
                % metric 1: 
                %   contribution based on ANOVA
                anovaTable = anova(model_4var_amide);
                contri_amide_all(:,idxalgo+(idxprepro-1)*3) = 100 * anovaTable.SumSq(1:4) / sum(anovaTable.SumSq); 
            
                anovaTable = anova(model_4var_guan);
                contri_guan_all(:,idxalgo+(idxprepro-1)*3) = 100 * anovaTable.SumSq(1:4) / sum(anovaTable.SumSq);
            
                anovaTable = anova(model_4var_noe);
                contri_noe_all(:,idxalgo+(idxprepro-1)*3) = 100 * anovaTable.SumSq(1:4) / sum(anovaTable.SumSq);
            
                % metric 2&3: 
                %   sensitivity based on linear regression slope, unit [%/mM]
                %   overall R2
                DeltaZ_guan = model_4var_guan.Variables{:,end};
                DeltaZ_amide = model_4var_amide.Variables{:,end};
                DeltaZ_noe = model_4var_noe.Variables{:,end};
            
                linFit = fitlm(model_4var_guan.Variables{:,1},DeltaZ_guan);
                sensR2_guan_all(1,idxalgo+(idxprepro-1)*3) = linFit.Coefficients.Estimate(2); % [%/mM]
        
                linFit = fitlm(model_4var_amide.Variables{:,2},DeltaZ_amide);
                sensR2_amide_all(1,idxalgo+(idxprepro-1)*3) = linFit.Coefficients.Estimate(2);
        
                linFit = fitlm(model_4var_noe.Variables{:,4},DeltaZ_noe);
                sensR2_noe_all(1,idxalgo+(idxprepro-1)*3) = linFit.Coefficients.Estimate(2);

                % metric 3:
                %   overall R2
                if idxalgo == 2 || idxalgo == 5 % LDA
                    R2flag = 0; % average positive and negative respectively
                else
                    R2flag = 1;
                end
                sensR2_amide_all(2,idxalgo+(idxprepro-1)*3) = calR2(model_4var_amide.Variables{:,2}, model_4var_amide.Variables{:,end},R2flag);
                sensR2_guan_all(2,idxalgo+(idxprepro-1)*3) = calR2(model_4var_guan.Variables{:,1}, model_4var_guan.Variables{:,end},R2flag);
                sensR2_noe_all(2,idxalgo+(idxprepro-1)*3) = calR2(model_4var_noe.Variables{:,4}, model_4var_noe.Variables{:,end},R2flag);

                % metric 4: difference with steady-state (Ts->infty)
                if idxTs == 1
                    % 1st Ts for steady-state (Ts->infty)
                    ssContrast_all(1,idxalgo+(idxprepro-1)*3,:) = Zamide_DZ_vec(:);
                    ssContrast_all(2,idxalgo+(idxprepro-1)*3,:) = Zguan_DZ_vec(:);
                    ssContrast_all(3,idxalgo+(idxprepro-1)*3,:) = Znoe_DZ_vec(:);
                end
                diffss_meanList_all_givenB1(1,idxalgo+(idxprepro-1)*3,idxTs) = ...
                    mean(ssContrast_all(1,idxalgo+(idxprepro-1)*3,:) - Zamide_DZ_vec(:),'all');
                diffss_meanList_all_givenB1(2,idxalgo+(idxprepro-1)*3,idxTs) = ...
                    mean(ssContrast_all(2,idxalgo+(idxprepro-1)*3,:) - Zguan_DZ_vec(:),'all');
                diffss_meanList_all_givenB1(3,idxalgo+(idxprepro-1)*3,idxTs) = ...
                    mean(ssContrast_all(3,idxalgo+(idxprepro-1)*3,:) - Znoe_DZ_vec(:),'all');


            end
        end

        zspec_vec_givenB1(:,:,idxTs) = zspec_vec;
        zspecQUASS_vec_givenB1(:,:,idxTs) = zspecQUASS_vec;
        contri_amide_all_givenB1(:,:,idxTs) = contri_amide_all;
        contri_guan_all_givenB1(:,:,idxTs) = contri_guan_all;
        contri_noe_all_givenB1(:,:,idxTs) = contri_noe_all;
        sensR2_amide_all_givenB1(:,:,idxTs) = sensR2_amide_all;
        sensR2_guan_all_givenB1(:,:,idxTs) = sensR2_guan_all;
        sensR2_noe_all_givenB1(:,:,idxTs) = sensR2_noe_all;

        fprintf(repmat('\b',1,backNum));
    end

    %% save data
    simData_B1Path = fullfile(outputfolder, "simData_B1_"+num2str(B1List(idxB1))+".mat");
    save(simData_B1Path,...
        'offs_trunc','zspec_vec_givenB1','zspecQUASS_vec_givenB1',...
        'contri_amide_all_givenB1','sensR2_amide_all_givenB1',...
        'contri_guan_all_givenB1','sensR2_guan_all_givenB1',...
        'contri_noe_all_givenB1','sensR2_noe_all_givenB1',...
        'diffss_meanList_all_givenB1');
    fprintf(" finished on "+string(datetime('now','InputFormat','yyyy-mm-dd HH:MM:SS'))+"\n");

end
delete(gcp('nocreate'));

function R2 = calR2(X, Y, R2flag)
    warning off
    if ~isempty(find(Y<0, 1)) && R2flag == 0
        % for LDA, the reference signal would be influenced by MT effect,
        % which is prominent when B1 is high, resulting negative Lorentzian
        % difference. In such case, a seperate R2 calculation process is
        % employed.
        Yposidx = find(Y>=0);
        Ynegidx = find(Y<0);
    
        linFit_pos = fitlm(X(Yposidx),Y(Yposidx));
        linFit_neg = fitlm(X(Ynegidx),Y(Ynegidx));
    
        R2_pos = linFit_pos.Rsquared.Ordinary;
        R2_neg = linFit_neg.Rsquared.Ordinary;

        if isnan(R2_pos)
            R2_pos = 0;
        end
        if isnan(R2_neg)
            R2_neg = 0;
        end
        R2 = R2_pos * length(Yposidx)/length(Y) + R2_neg * length(Ynegidx)/length(Y);
    else
        linFit = fitlm(X, Y);
        R2 = linFit.Rsquared.Ordinary;
    end
end
