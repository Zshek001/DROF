% perform PLOF on phantom data (random setting)
% use multi-linear regression to assess dependence on MTf, NOEf
% 2025-0121
clear
close all
warning off

addpath ..\toolbox\
addpath ..\toolbox\fitting_main\

display_on = 1;
multistart_N = 5; % times for MultiStart
fileName = ".\simData_random5000.mat";
outputfolder = ".\ran5000-fig3\";

if ~exist("outputfolder",'dir')
    mkdir(outputfolder)
end

%% load data
load(fileName,"offs","zspec"); % m0 not included
    % offs: [nf,1]
    %   offs(1): m0
    % zspec: [nf, Npixel]
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

%% (1) Delta Z
processing_DZ_methods = {@fitting_MTR,@fitting_LDA,@fitting_LDA_2pool,@fitting_MPLF,@fitting_MPLF_2step};
output_filenames = {'DZ_MTR','DZ_LDA','DZ_LDA_2pool','CO_MPLF','CO_MPLF_2step'};
for i = 1:length(processing_DZ_methods)
    method = processing_DZ_methods{i};
    [Zamide_DZ_vec, Zguan_DZ_vec] = method(offs, zspec_vec, multistart_N);

    % construct response variable vector, [nPixel,1]
    responVar_guan = Zguan_DZ_vec(:);
    responVar_amide = Zamide_DZ_vec(:);

    % linear regression
    model_4var_guan = fitlm(predictVar, responVar_guan,'ResponseVar','metric_guan','PredictorVars',{'guanf','amidef','MTf','NOEf'});
    model_4var_amide = fitlm(predictVar, responVar_amide,'ResponseVar','metric_amide','PredictorVars',{'guanf','amidef','MTf','NOEf'});


    save(outputfolder+"output_"+output_filenames{i}+".mat",'model_4var_amide','model_4var_guan');

    fprintf('saved to %s\n', outputfolder+"output_"+output_filenames{i}+".mat");
end

%% (2) R1rho
B0 = 3; % T
FitParam.R1 = 1; % brain ~ 1s @ 3T
FitParam.satpwr = 0.8; % [uT]
FitParam.tsat = 2; % [s]
FitParam.Magfield = 42.5764 * B0; % [Hz]
FitParam.PeakOffset = 3.5; % used in PLOF, position to be polynomialize

processing_R1rho_methods = {@fitting_SROF,@fitting_DROF,@fitting_PLOF};
output_filenames = {'R1rho_SROF','R1rho_DROF','R1rho_PLOF'};
for i = 1:length(processing_R1rho_methods)
    method = processing_R1rho_methods{i};

    nOut = nargout(method);
    if nOut == 3
        [~,~, fit_para] = method(offs, zspec_vec, FitParam, multistart_N); 
    elseif nOut ==4
        [~,~,~, fit_para] = method(offs, zspec_vec, FitParam, multistart_N); 
    end

    Zamide_DZ_vec = zeros(Npixel,1);
    Zguan_DZ_vec = zeros(Npixel,1);
    if strcmp(output_filenames{i},'R1rho_PLOF')
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
    
    save(outputfolder+"output_"+output_filenames{i}+".mat",'model_4var_amide','model_4var_guan','model_4var_amide_DZ','model_4var_guan_DZ','Zamide_DZ_vec','Zguan_DZ_vec','fit_para');

    fprintf('saved to %s\n', outputfolder+"output_"+output_filenames{i}+".mat");
end










