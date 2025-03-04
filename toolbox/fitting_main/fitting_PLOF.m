function [Zamide_DZ_vec, Zguan_DZ_vec, fit_para] = fitting_PLOF(offs, zspec_vec, FitParam, multistart_N)
% OUTPUT
%   fit_para: [8, Npixel]
%       Water (1-2), Amide (3-5), Guan (6-8)
    %% PLOF setting
    % saturation settings
    % pulse1_pwr = 0.8; % in uT
    % pulse1_dur = 2; % pulse duration in s
    % gamma_hz = 42.5764;
    % B0 = 3;
    
    % FitParam.CalSNR = 0; % use m0 img to calculate SNR
    % FitParam.SNRimg = zspec_vec(1,:); 
    FitParam.PeakOffset = 3.5; % used in CurveFunction, position to be polynomialize
    % FitParam.R1 = 1; % no T1 correction
    % FitParam.satpwr = pulse1_pwr; % saturation power (uT)
    % FitParam.tsat = pulse1_dur; % saturation length (second), 100s to make it steady-state
    % FitParam.Magfield = gamma_hz * B0;
    FitParam.MultiStart = multistart_N;
    WholeRange = [0.5, 8];    % signal for background/CEST fitting
    PeakRange = [1, 6];     % extract background signal
    
    [nf,Npixel] = size(zspec_vec);
    Ramide_vec = zeros(1,Npixel);
    Rguan_vec = zeros(1,Npixel);
    
    Zamide_DZ_vec = zeros(1,Npixel);
    Zguan_DZ_vec = zeros(1,Npixel);
    Zamide_Co_vec = zeros(1,Npixel);
    Zguan_Co_vec = zeros(1,Npixel);
    fit_para = []; % [npool, Npixel]
    
    MTbg_vec = zeros(1,Npixel);
    
    %% PLOF fitting
    backNum = 0;
    tic
    for i = 1:Npixel
        fprintf(repmat('\b',1,backNum));
        backNum = fprintf('Calculation process: %d/%d',i,Npixel);
                  
        Z_spectrum = zspec_vec(:,i);
        Z_spectrum = squeeze(Z_spectrum);
    
        % only use zspec without m0
        [FitResult,FitParam] = PLOF_custom(offs, Z_spectrum, WholeRange, PeakRange, FitParam);
    
        Ramide_vec(i) = 1000*FitResult.DeltaRpeak1;
        Rguan_vec(i) = 1000*FitResult.DeltaRpeak2;
        MTbg_vec(i) = 100 * (1 - FitResult.MTbackground);
    
        Zamide_DZ_vec(i) = 100*FitResult.DeltaZpeak1; % use delta Z as metric
        Zguan_DZ_vec(i) = 100*FitResult.DeltaZpeak2;
    
        Zamide_Co_vec(i) = 100*FitResult.Rpeak1; % use Lorentzian coeff as metric
        Zguan_Co_vec(i) = 100*FitResult.Rpeak2;

        fit_para = [fit_para,reshape(FitResult.Coefficents,[],1)];

    end
    fprintf(repmat('\b',1,backNum));
    fprintf('elapsed time: %.4f s\n', toc);
end