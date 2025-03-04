function [FitResult_best,FitResult_wor, FitParam] = PLOF_custom_beworst(Offset, Saturation, WholeRange, PeakRange, FitParam)
% Fit the Z-spectrum using polynomial and Lorentzian line-shape fitting (PLOF) method
% Please contact Kexin Wang at kwang101@jhu.edu if you have any questions about the code. 
% modified from 2024-0923

warning off

Offset = reshape(Offset,[],1);
Saturation = reshape(Saturation,[],1);

% Z-spectra ROI
[~, idx1] = min(abs(Offset - min(WholeRange)));
[~, idx2] = min(abs(Offset - max(WholeRange)));
idxlow = min([idx1,idx2]); idxupp = max([idx1,idx2]);
Offset = Offset(idxlow:idxupp);
Saturation = Saturation(idxlow:idxupp);

[~, idx1] = min(abs(Offset - min(PeakRange)));
[~, idx2] = min(abs(Offset - max(PeakRange)));
idxlow = min([idx1,idx2]); idxupp = max([idx1,idx2]);
Offset_background = Offset; Offset_background(idxlow:idxupp) = [];
Saturation_background = Saturation; Saturation_background(idxlow:idxupp) = [];

FitResult_best.Offset = Offset;
FitResult_best.Saturation = Saturation;
FitResult_best.Offset_background = Offset_background;
FitResult_best.Saturation_background = Saturation_background;

%-----------------------------------------------%
%-------------- Two-step fitting ---------------%
%-----------------------------------------------%

FitResult_best.xindex = reshape(linspace(min(Offset), max(Offset), 100),[],1); % interpolated offs, in ppm

% (1) first-order polynomial
%               1. Background                  2. Amide                   3. Guan
%       A1      G1      P0      P1        A2      G2      dw2         A3      G3      dw3
x0 = [  5       1       0.5     -19       0.1     1       3.5         0.1     1       2];
lb = [  0       0       0       -100      1e-8    0.001   3           1e-8    0.01    1];
ub = [  20      100     100     0         20      2       4           10      2       3];

% % (2) second-order polynomial
% %                   1. Background                      2. Amide                   3. Guan
%       A1      G1      P0      P1      P2        A2      G2      dw2         A3      G3      dw3
% x0 = [  0.1     1       0.5     -19     1         0.1     1       3.5         0.1     1       2];
% lb = [  0       0       0       -100    -10       1e-8    0.001   3           1e-8    0.01    1];
% ub = [  10      100     100     0       10        10      2       4           10      2       3];
    
%-------------step 1: background fitting --------------%

FitParam.peak = 0;
%   peak = 0: only fit background;
%        = 1: all peaks are used for the fitting;
%        = 2: only backgound and peak_1 are used for the fitting;
%        = 3: only backgound and peak_2 are used for the fitting;

options = optimset('MaxFunEvals', 1e6,'TolFun', 1e-6,'TolX', 1e-6, 'Display',  'off' );
CurveFunc_2var = @(x,xdata) CurveFunction(x,xdata,FitParam);
% [FitResult.Coefficents, resnorm] = lsqcurvefit(CurveFunc_2var, x0, Offset_background, Saturation_background, lb, ub, options);
problem = createOptimProblem('lsqcurvefit','x0',x0,'objective',CurveFunc_2var,...
    'lb',lb,'ub',ub,'xdata',Offset_background,'ydata',Saturation_background);

rng default % for reproducibility
ms = MultiStart('Display','off','OutputFcn', @record_max_fval_params);
[FitResult_best.Coefficents, ~] = run(ms,problem,FitParam.MultiStart);
FitResult_wor.Coefficents = record_max_fval_params();

% fitting result of Z-spectrum and R-spectrum
[FitResult_best.Background, Rbak] = CurveFunction(FitResult_best.Coefficents, FitResult_best.xindex, FitParam);
Backgroundtmp = CurveFunction(FitResult_best.Coefficents, Offset, FitParam);
FitResult_best.DeltaZspectrum = Backgroundtmp - Saturation; 

% R^2 for the evaluation of fitting
FitResult_best.RsquareBG = 1 - goodnessOfFit(CurveFunction(FitResult_best.Coefficents, Offset_background, FitParam), Saturation_background, 'NMSE');

%   fitted background at 3.5 ppm
[res, idx] = min(abs(FitResult_best.xindex - 3.5)); % will not influence result
FitResult_best.MTbackground = FitResult_best.Background(idx); 

%-------------step 2: CEST all-peak fitting --------------%
switch length(x0)
    case 10 % first-order
        ind_poly = 1:4;
    case 11 % second-order
        ind_poly = 1:5;
end
        

FitParam.peak = 1;
%   peak = 0: only fit background;
%        = 1: all peaks are used for the fitting;
%        = 2: only backgound and peak_1 are used for the fitting;
%        = 3: only backgound and peak_2 are used for the fitting;

%% worst
% fix the background parameters
lb(ind_poly) = FitResult_wor.Coefficents(ind_poly); 
ub(ind_poly) = FitResult_wor.Coefficents(ind_poly); 

% [FitResult.Coefficents, resnorm] = lsqcurvefit(@CurveFunction, x0, Offset, Saturation, lb, ub, options, FitParam);
CurveFunc_2var = @(x,xdata) CurveFunction(x,xdata,FitParam);
problem = createOptimProblem('lsqcurvefit','x0',x0,'objective',CurveFunc_2var,...
    'lb',lb,'ub',ub,'xdata',Offset,'ydata',Saturation);
ms = MultiStart('Display','off','OutputFcn', @record_max_fval_params);
[FitResult_wor.Coefficents, resnorm] = run(ms,problem,FitParam.MultiStart);
FitResult_wor.Coefficents = record_max_fval_params();

%% best
% fix the background parameters
lb(ind_poly) = FitResult_best.Coefficents(ind_poly); 
ub(ind_poly) = FitResult_best.Coefficents(ind_poly); 

% [FitResult.Coefficents, resnorm] = lsqcurvefit(@CurveFunction, x0, Offset, Saturation, lb, ub, options, FitParam);
CurveFunc_2var = @(x,xdata) CurveFunction(x,xdata,FitParam);
problem = createOptimProblem('lsqcurvefit','x0',x0,'objective',CurveFunc_2var,...
    'lb',lb,'ub',ub,'xdata',Offset,'ydata',Saturation);
ms = MultiStart('Display','off','OutputFcn', @record_max_fval_params);
[FitResult_best.Coefficents, resnorm] = run(ms,problem,FitParam.MultiStart);

[FitResult_best.Curve, Rtmp] = CurveFunction(FitResult_best.Coefficents, FitResult_best.xindex, FitParam);

FitResult_best.DeletaFitZ = FitResult_best.Background-FitResult_best.Curve;
FitResult_best.DeletaFitR = Rtmp - Rbak;


%------------ Calculate deltaZ and assign the Coefficents ---------------%

FitResult_best.Rpeak1 = FitResult_best.Coefficents(length(ind_poly)+1);
FitResult_best.FitPeakOffset1 = FitResult_best.Coefficents(length(ind_poly)+3);
FitResult_best.Rpeak2 = FitResult_best.Coefficents(length(ind_poly)+4);
FitResult_best.FitPeakOffset2 = FitResult_best.Coefficents(length(ind_poly)+6); 

FitResult_best.MT = FitResult_best.Coefficents(1); 
FitResult_best.RsquareAll = 1 - goodnessOfFit(CurveFunction(FitResult_best.Coefficents, Offset, FitParam), Saturation, 'NMSE');


%------------ Calculate DeltaZ-spectrum ---------------%
flag = 1;
% amide peak
FitParam.peak = 2;

[Ztmp, Rtmp]= CurveFunction(FitResult_best.Coefficents, FitResult_best.xindex, FitParam); 
FitResult_best.ZamideFit = FitResult_best.Background -Ztmp;
FitResult_best.RamideFit = Rtmp - Rbak;
switch flag
    case 1
        % (1) maximum peak
        [FitResult_best.DeltaZpeak1, idx] = max(FitResult_best.ZamideFit); 
        FitResult_best.DeltaZpeak1Offset = FitResult_best.xindex(idx); 

        [FitResult_best.DeltaRpeak1, idx] = max(FitResult_best.RamideFit);
        FitResult_best.DeltaRpeak1Offset = FitResult_best.xindex(idx); 
        
    case 2
        % (2) fixed peak
        offs_amide = 3.5;
        [~, idx] = min(abs(Offset - offs_amide));
        FitResult_best.DeltaZpeak1 = FitResult_best.ZamideFit(idx);
        FitResult_best.DeltaRpeak1 = FitResult_best.RamideFit(idx);

        FitResult_best.DeltaZpeak1Offset = FitResult_best.xindex(idx); 
        FitResult_best.DeltaRpeak1Offset = FitResult_best.xindex(idx); 
end

% guan peak
FitParam.peak = 3;

[Ztmp, Rtmp] = CurveFunction(FitResult_best.Coefficents, FitResult_best.xindex, FitParam); 
FitResult_best.ZguanFit =FitResult_best.Background -Ztmp;
FitResult_best.RguanFit =Rtmp - Rbak;
switch flag
    case 1
        % (1) maximum peak
        [FitResult_best.DeltaZpeak2, idx] = max(FitResult_best.ZguanFit); 
        FitResult_best.DeltaZpeak2Offset = FitResult_best.xindex(idx); 

        [FitResult_best.DeltaRpeak2, idx] = max(FitResult_best.RguanFit);
        FitResult_best.DeltaRpeak2Offset = FitResult_best.xindex(idx); 
        
    case 2
        % (2) fixed peak
        offs_guan = 2;
        [~, idx] = min(abs(Offset - offs_guan));
        FitResult_best.DeltaZpeak2 = FitResult_best.ZguanFit(idx);
        FitResult_best.DeltaRpeak2 = FitResult_best.RguanFit(idx);

        FitResult_best.DeltaZpeak2Offset = FitResult_best.xindex(idx); 
        FitResult_best.DeltaRpeak2Offset = FitResult_best.xindex(idx); 
end

end

function stop = record_max_fval_params(optimValues, state)
    persistent max_fval_store worst_params_store;
    stop = false;

    if nargin == 0
        % max_fval = max_fval_store;
        stop = worst_params_store;
        return;
    end
    
    if strcmp(state, 'init')
        max_fval_store = -Inf;
        worst_params_store = [];
    elseif strcmp(state, 'iter')
        if optimValues.localsolution.Fval > max_fval_store
            max_fval_store = optimValues.localsolution.Fval;
            worst_params_store = optimValues.localsolution.X;
        end
    end
end