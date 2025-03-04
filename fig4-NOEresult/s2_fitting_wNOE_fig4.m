% perform PLOF on phantom data (random setting)
% use multi-linear regression to assess dependence on MTf, NOEf
%       compare amide, guan, NOE
%       compare two offs
%       compare MTR, LDA, MPLF, MPLF2-ss
% 2025-0121
%   use user19.txt offset (55p)
%   save all delta Z
clear
% close all
warning off

addpath ..\toolbox\
addpath ..\toolbox\fitting_main\

display_on = 1;
multistart_N = 5; % times for MultiStart
fileName = ".\simData_random5000.mat";
outputfolder = ".\ran5000-fig4\";

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
processing_DZ_methods = {@fitting_MTR,@fitting_LDA,@fitting_MPLF};
output_filenames = {'DZ_MTR','DZ_LDA','CO_MPLF'};
for i = 1:length(processing_DZ_methods)
    method = processing_DZ_methods{i};
    [Zamide_DZ_vec, Zguan_DZ_vec, Znoe_DZ_vec] = method(offs, zspec_vec, multistart_N);

    % construct response variable vector, [nPixel,1]
    responVar_guan = Zguan_DZ_vec(:);
    responVar_amide = Zamide_DZ_vec(:);
    responVar_noe = Znoe_DZ_vec(:);

    % linear regression
    model_4var_guan = fitlm(predictVar, responVar_guan,'ResponseVar','metric_guan','PredictorVars',{'guanf','amidef','MTf','NOEf'});
    model_4var_amide = fitlm(predictVar, responVar_amide,'ResponseVar','metric_amide','PredictorVars',{'guanf','amidef','MTf','NOEf'});
    model_4var_noe = fitlm(predictVar, responVar_noe,'ResponseVar','metric_amide','PredictorVars',{'guanf','amidef','MTf','NOEf'});

    save(outputfolder+"output_"+output_filenames{i}+".mat",'model_4var_amide','model_4var_guan','model_4var_noe');

    fprintf('saved to %s\n', outputfolder+"output_"+output_filenames{i}+".mat");
end

%% (2) R1rho
B0 = 3; % T
FitParam.R1 = 1; % brain ~ 1s @ 3T
FitParam.satpwr = 0.8; % [uT]
FitParam.tsat = 2; % [s]
FitParam.Magfield = 42.5764 * B0; % [Hz]
FitParam.PeakOffset = 3.5; % used in PLOF, position to be polynomialize

processing_R1rho_methods = {@fitting_DROF};
output_filenames = {'R1rho_DROF'};
for i = 1:length(processing_R1rho_methods)
    method = processing_R1rho_methods{i};
    [~,~,~, fit_para] = method(offs, zspec_vec, FitParam, multistart_N); 

    Zamide_DZ_vec = zeros(Npixel,1);
    Zguan_DZ_vec = zeros(Npixel,1);
    
    indamide = 5; indguan = 14; indnoe = 8;
    % calculate zspec for analysis
    R1rho_ss = @(par,offs) lorentzMultipool_R1rho(par,offs,FitParam);
    for i_pixel = 1:Npixel
        par = fit_para(:,i_pixel);
        Z_dsMT = R1rho_ss([par([1:4,11:13])],offs); % Water + MT
        Z_bak = R1rho_ss([par([1:4,8:13])],offs); % Water + MT + NOE
        Z_guan = R1rho_ss(par([1:4,8:16]), offs); % Water + MT + NOE + guan
        Z_amide = R1rho_ss(par(1:13), offs); % Water + MT + NOE + amide
        Z_ful = R1rho_ss(par, offs); % Water + MT + NOE + guan + amide
        Zamide_DZ_vec(i_pixel) = 100*max(Z_bak-Z_amide);
        Zguan_DZ_vec(i_pixel) = 100*max(Z_bak-Z_guan);
        Znoe_DZ_vec(i_pixel) = 100*max(Z_dsMT-Z_bak);
    end

    %% save
    Zamide_Co_vec = 1000*fit_para(indamide,:); % unit ms^-1
    Zguan_Co_vec = 1000*fit_para(indguan,:);
    Znoe_Co_vec = 1000*fit_para(indnoe,:);

    % construct response variable vector, [nPixel,1]
    responVar_guan = Zguan_Co_vec(:);
    responVar_amide = Zamide_Co_vec(:);
    responVar_noe = Znoe_Co_vec(:);
    
    % linear regression
    model_4var_guan = fitlm(predictVar, responVar_guan,'ResponseVar','metric_guan','PredictorVars',{'guanf','amidef','MTf','NOEf'});
    model_4var_amide = fitlm(predictVar, responVar_amide,'ResponseVar','metric_amide','PredictorVars',{'guanf','amidef','MTf','NOEf'});
    model_4var_noe = fitlm(predictVar, responVar_noe,'ResponseVar','metric_amide','PredictorVars',{'guanf','amidef','MTf','NOEf'});

    responVar_guan = Zguan_DZ_vec(:);
    responVar_amide = Zamide_DZ_vec(:);
    responVar_noe = Znoe_DZ_vec(:);
    model_4var_guan_DZ = fitlm(predictVar, responVar_guan,'ResponseVar','metric_guan','PredictorVars',{'guanf','amidef','MTf','NOEf'});
    model_4var_amide_DZ = fitlm(predictVar, responVar_amide,'ResponseVar','metric_amide','PredictorVars',{'guanf','amidef','MTf','NOEf'});
    model_4var_noe_DZ = fitlm(predictVar, responVar_noe,'ResponseVar','metric_amide','PredictorVars',{'guanf','amidef','MTf','NOEf'});


    save(outputfolder+"output_"+output_filenames{i}+".mat",'model_4var_amide','model_4var_guan','model_4var_noe','model_4var_amide_DZ','model_4var_guan_DZ','model_4var_noe_DZ','Zamide_DZ_vec','Zguan_DZ_vec','Znoe_DZ_vec','fit_para');


    fprintf('saved to %s\n', outputfolder+"output_"+output_filenames{i}+".mat");
end
