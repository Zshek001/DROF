% use MultiStart(50), store best and worst fit_para
%       for MPLF, MPLF2, SROF, DROF, PLOF
% 2025-0121 
%   save all delta Z
%   use random dataset
clear
close all
warning off

addpath ..\toolbox\
addpath ..\toolbox\fitting_main\

display_on = 1;
multistart_N = 50; % times for MultiStart
fileName = ".\simData_random5000.mat";
outputfolder = ".\BesAWor-ran5000-MS50\";

if ~exist("outputfolder",'dir')
    mkdir(outputfolder)
end

%% load data
load(fileName,"offs","zspec"); % m0 not included
    % offs: [nf,1]
    %   offs(1): m0
    % zspec: [nf, nguanf, nguank, namidef]
load(fileName,"paraList",'paraInfo');
    % paraList: [7,Npixel]
    %   MTf_relative, Guank_Hz, Amidef_mM, Guanf_mM, NOEf_mM, Amidef_Hz, NOEk_Hz

MTf_mM = paraList(1,:);
% guank = paraList(2,:); % Hz
amidef_mM = paraList(3,:); % mM
guanf_mM = paraList(4,:);
noef_mM = paraList(5,:);

[nf,Npixel] = size(zspec);
zspec_vec = zspec;

% construct predictor variable matrix for linear regression, [nPixel,nvar]
predictVar = [guanf_mM(:),amidef_mM(:),MTf_mM(:),noef_mM(:)];

%% (1) Delta Z best
processing_DZ_methods = {@fitting_MPLF_beworst,@fitting_MPLF_2step_beworst};
output_filenames = {'CO_MPLF','CO_MPLF_2step'};
for i = 1:length(processing_DZ_methods)
    method = processing_DZ_methods{i};
    [Zamide_DZ_best, Zguan_DZ_best, Zamide_DZ_wor, Zguan_DZ_wor] = method(offs, zspec_vec, multistart_N);

    %% best
    % construct response variable vector, [nPixel,1]
    responVar_guan = Zguan_DZ_best(:);
    responVar_amide = Zamide_DZ_best(:);

    % linear regression
    model_4var_guan = fitlm(predictVar, responVar_guan,'ResponseVar','metric_guan','PredictorVars',{'guanf','amidef','MTf','NOEf'});
    model_4var_amide = fitlm(predictVar, responVar_amide,'ResponseVar','metric_amide','PredictorVars',{'guanf','amidef','MTf','NOEf'});

    save(outputfolder+"output_"+output_filenames{i}+"_best.mat",'model_4var_amide','model_4var_guan');

    %% worst
    % construct response variable vector, [nPixel,1]
    responVar_guan = Zguan_DZ_wor(:);
    responVar_amide = Zamide_DZ_wor(:);

    % linear regression
    model_4var_guan = fitlm(predictVar, responVar_guan,'ResponseVar','metric_guan','PredictorVars',{'guanf','amidef','MTf','NOEf'});
    model_4var_amide = fitlm(predictVar, responVar_amide,'ResponseVar','metric_amide','PredictorVars',{'guanf','amidef','MTf','NOEf'});

    save(outputfolder+"output_"+output_filenames{i}+"_worst.mat",'model_4var_amide','model_4var_guan');

    fprintf('saved to %s\n', outputfolder+"output_"+output_filenames{i}+".mat");
end

%% (2) R1rho best
B0 = 3; % T
FitParam.R1 = 1; % brain ~ 1s @ 3T
FitParam.satpwr = 0.8; % [uT]
FitParam.tsat = 2; % [s]
FitParam.Magfield = 42.5764 * B0; % [Hz]
FitParam.PeakOffset = 3.5; % used in PLOF, position to be polynomialize

processing_R1rho_methods = {@fitting_SROF_beworst,@fitting_DROF_beworst,@fitting_PLOF_beworst};
output_filenames = {'DZ_SROF','DZ_DROF','DZ_PLOF'};
for i = 1:length(processing_R1rho_methods)
    method = processing_R1rho_methods{i};
    [~,~, fit_para_best,fit_para_wor] = method(offs, zspec_vec, FitParam, multistart_N); 

    %% best
    fit_para = fit_para_best;
    Zamide_DZ_vec = zeros(Npixel,1);
    Zguan_DZ_vec = zeros(Npixel,1);
    if strcmp(output_filenames{i},'DZ_PLOF')
        nterm = size(fit_para,1);
        indamide = nterm - 5;
        indguan = nterm - 2;
        % calculate zspec for analysis
        WholeRange = [0.5,8];
        [~, idx1] = min(abs(offs - min(WholeRange)));
        [~, idx2] = min(abs(offs - max(WholeRange)));
        idxlow = min([idx1,idx2]); idxupp = max([idx1,idx2]);
        offs_PLOF = offs(idxlow:idxupp);
        for i_pixel = 1:Npixel
            par = fit_para(:,i_pixel);
            FitParam.peak = 0; Z_bak = CurveFunction(par,offs_PLOF,FitParam);
            FitParam.peak = 1; Z_ful = CurveFunction(par,offs_PLOF,FitParam);
            FitParam.peak = 2; Z_amide = CurveFunction(par,offs_PLOF,FitParam);
            FitParam.peak = 3; Z_guan = CurveFunction(par,offs_PLOF,FitParam);
            Zamide_DZ_vec(i_pixel) = 100*max(Z_bak-Z_amide);
            Zguan_DZ_vec(i_pixel) = 100*max(Z_bak-Z_guan);
        end
    else
        indamide = 5; indguan = 14;
        % calculate zspec for analysis
        R1rho_ss = @(par,offs) lorentzMultipool_R1rho(par,offs,FitParam);
        for i_pixel = 1:Npixel
            par = fit_para(:,i_pixel);
            Z_bak = R1rho_ss([par([1:4,8:13])],offs); % Water + MT + NOE
            Z_guan = R1rho_ss(par([1:4,8:16]), offs); % Water + MT + NOE + guan
            Z_amide = R1rho_ss(par(1:13), offs); % Water + MT + NOE + amide
            Z_ful = R1rho_ss(par, offs); % Water + MT + NOE + guan + amide
            Zamide_DZ_vec(i_pixel) = 100*max(Z_bak-Z_amide);
            Zguan_DZ_vec(i_pixel) = 100*max(Z_bak-Z_guan);
        end
    end
    Zamide_Co_vec = 1000*fit_para(indamide,:); % unit ms^-1
    Zguan_Co_vec = 1000*fit_para(indguan,:);

    % construct response variable vector, [nPixel,1]
    responVar_guan = Zguan_Co_vec(:);
    responVar_amide = Zamide_Co_vec(:);
    
    % linear regression
    model_4var_guan = fitlm(predictVar, responVar_guan,'ResponseVar','metric_guan','PredictorVars',{'guanf','amidef','MTf','NOEf'});
    model_4var_amide = fitlm(predictVar, responVar_amide,'ResponseVar','metric_amide','PredictorVars',{'guanf','amidef','MTf','NOEf'});

    % DZ
    responVar_guan = Zguan_DZ_vec(:);
    responVar_amide = Zamide_DZ_vec(:);
    model_4var_guan_DZ = fitlm(predictVar, responVar_guan,'ResponseVar','metric_guan','PredictorVars',{'guanf','amidef','MTf','NOEf'});
    model_4var_amide_DZ = fitlm(predictVar, responVar_amide,'ResponseVar','metric_amide','PredictorVars',{'guanf','amidef','MTf','NOEf'});
    
    save(outputfolder+"output_"+output_filenames{i}+"_best.mat",'model_4var_amide','model_4var_guan','model_4var_amide_DZ','model_4var_guan_DZ','fit_para');
    
    %% worst
    fit_para = fit_para_wor;
    Zamide_DZ_vec = zeros(Npixel,1);
    Zguan_DZ_vec = zeros(Npixel,1);
    if strcmp(output_filenames{i},'DZ_PLOF')
        nterm = size(fit_para,1);
        indamide = nterm - 5;
        indguan = nterm - 2;
        % calculate zspec for analysis
        WholeRange = [0.5,8];
        [~, idx1] = min(abs(offs - min(WholeRange)));
        [~, idx2] = min(abs(offs - max(WholeRange)));
        idxlow = min([idx1,idx2]); idxupp = max([idx1,idx2]);
        offs_PLOF = offs(idxlow:idxupp);
        for i_pixel = 1:Npixel
            par = fit_para(:,i_pixel);
            FitParam.peak = 0; Z_bak = CurveFunction(par,offs_PLOF,FitParam);
            FitParam.peak = 1; Z_ful = CurveFunction(par,offs_PLOF,FitParam);
            FitParam.peak = 2; Z_amide = CurveFunction(par,offs_PLOF,FitParam);
            FitParam.peak = 3; Z_guan = CurveFunction(par,offs_PLOF,FitParam);
            Zamide_DZ_vec(i_pixel) = 100*max(Z_bak-Z_amide);
            Zguan_DZ_vec(i_pixel) = 100*max(Z_bak-Z_guan);
        end
    else
        indamide = 5; indguan = 14;
        % calculate zspec for analysis
        R1rho_ss = @(par,offs) lorentzMultipool_R1rho(par,offs,FitParam);
        for i_pixel = 1:Npixel
            par = fit_para(:,i_pixel);
            Z_bak = R1rho_ss([par([1:4,8:13])],offs); % Water + MT + NOE
            Z_guan = R1rho_ss(par([1:4,8:16]), offs); % Water + MT + NOE + guan
            Z_amide = R1rho_ss(par(1:13), offs); % Water + MT + NOE + amide
            Z_ful = R1rho_ss(par, offs); % Water + MT + NOE + guan + amide
            Zamide_DZ_vec(i_pixel) = 100*max(Z_bak-Z_amide);
            Zguan_DZ_vec(i_pixel) = 100*max(Z_bak-Z_guan);
        end
    end
    Zamide_Co_vec = 1000*fit_para(indamide,:); % unit ms^-1
    Zguan_Co_vec = 1000*fit_para(indguan,:);

    % construct response variable vector, [nPixel,1]
    responVar_guan = Zguan_Co_vec(:);
    responVar_amide = Zamide_Co_vec(:);
    
    % linear regression
    model_4var_guan = fitlm(predictVar, responVar_guan,'ResponseVar','metric_guan','PredictorVars',{'guanf','amidef','MTf','NOEf'});
    model_4var_amide = fitlm(predictVar, responVar_amide,'ResponseVar','metric_amide','PredictorVars',{'guanf','amidef','MTf','NOEf'});

    % DZ
    responVar_guan = Zguan_DZ_vec(:);
    responVar_amide = Zamide_DZ_vec(:);
    model_4var_guan_DZ = fitlm(predictVar, responVar_guan,'ResponseVar','metric_guan','PredictorVars',{'guanf','amidef','MTf','NOEf'});
    model_4var_amide_DZ = fitlm(predictVar, responVar_amide,'ResponseVar','metric_amide','PredictorVars',{'guanf','amidef','MTf','NOEf'});
    
    save(outputfolder+"output_"+output_filenames{i}+"_worst.mat",'model_4var_amide','model_4var_guan','model_4var_amide_DZ','model_4var_guan_DZ','fit_para');

    fprintf('saved to %s\n', outputfolder+"output_"+output_filenames{i}+".mat");
end
