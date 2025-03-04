function [Z_fit, Z_bak, Z_amide, Z_guan] = PLOF_custom_demo(Offset, Saturation, WholeRange, PeakRange, FitParam)
% Fit the Z-spectrum using polynomial and Lorentzian line-shape fitting (PLOF) method
% Please contact Kexin Wang at kwang101@jhu.edu if you have any questions about the code. 
% modified from 2024-0923
% 20241104

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

FitResult.Offset = Offset;
FitResult.Saturation = Saturation;
FitResult.Offset_background = Offset_background;
FitResult.Saturation_background = Saturation_background;

%-----------------------------------------------%
%-------------- Two-step fitting ---------------%
%-----------------------------------------------%

FitResult.xindex = reshape(linspace(min(Offset), max(Offset), 100),[],1); % interpolated offs, in ppm

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
ms = MultiStart('Display','off');
[FitResult.Coefficents, resnorm] = run(ms,problem,FitParam.MultiStart);

% fitting result of Z-spectrum and R-spectrum
[FitResult.Background, Rbak] = CurveFunction(FitResult.Coefficents, FitResult.xindex, FitParam);
Backgroundtmp = CurveFunction(FitResult.Coefficents, Offset, FitParam);
FitResult.DeltaZspectrum = Backgroundtmp - Saturation; 

% R^2 for the evaluation of fitting
FitResult.RsquareBG = 1 - goodnessOfFit(CurveFunction(FitResult.Coefficents, Offset_background, FitParam), Saturation_background, 'NMSE');

%   fitted background at 3.5 ppm
[res, idx] = min(abs(FitResult.xindex - 3.5)); % will not influence result
FitResult.MTbackground = FitResult.Background(idx); 

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

% fix the background parameters
lb(ind_poly) = FitResult.Coefficents(ind_poly); 
ub(ind_poly) = FitResult.Coefficents(ind_poly); 

% [FitResult.Coefficents, resnorm] = lsqcurvefit(@CurveFunction, x0, Offset, Saturation, lb, ub, options, FitParam);
CurveFunc_2var = @(x,xdata) CurveFunction(x,xdata,FitParam);
problem = createOptimProblem('lsqcurvefit','x0',x0,'objective',CurveFunc_2var,...
    'lb',lb,'ub',ub,'xdata',Offset,'ydata',Saturation);
ms = MultiStart('Display','off');
[FitResult.Coefficents, resnorm] = run(ms,problem,FitParam.MultiStart);

[FitResult.Curve, Rtmp] = CurveFunction(FitResult.Coefficents, FitResult.xindex, FitParam);

FitResult.DeletaFitZ = FitResult.Background-FitResult.Curve;
FitResult.DeletaFitR = Rtmp - Rbak;

%-------------20241104 curve as MPLF---------------------------------
FitParam.peak = 0;
Z_bak = CurveFunction(FitResult.Coefficents, Offset, FitParam);

FitParam.peak = 1;
Z_fit = CurveFunction(FitResult.Coefficents, Offset, FitParam);

FitParam.peak = 21;
Z_amide = CurveFunction(FitResult.Coefficents, Offset, FitParam);

FitParam.peak = 31;
Z_guan = CurveFunction(FitResult.Coefficents, Offset, FitParam);


%------------ Calculate deltaZ and assign the Coefficents ---------------%

FitResult.Rpeak1 = FitResult.Coefficents(length(ind_poly)+1);
FitResult.FitPeakOffset1 = FitResult.Coefficents(length(ind_poly)+3);
FitResult.Rpeak2 = FitResult.Coefficents(length(ind_poly)+4);
FitResult.FitPeakOffset2 = FitResult.Coefficents(length(ind_poly)+6); 

FitResult.MT = FitResult.Coefficents(1); 
FitResult.RsquareAll = 1 - goodnessOfFit(CurveFunction(FitResult.Coefficents, Offset, FitParam), Saturation, 'NMSE');


%------------ Calculate DeltaZ-spectrum ---------------%
flag = 1;
% amide peak
FitParam.peak = 2;

[Ztmp, Rtmp]= CurveFunction(FitResult.Coefficents, FitResult.xindex, FitParam); 
FitResult.ZamideFit = FitResult.Background -Ztmp;
FitResult.RamideFit = Rtmp - Rbak;
switch flag
    case 1
        % (1) maximum peak
        [FitResult.DeltaZpeak1, idx] = max(FitResult.ZamideFit); 
        FitResult.DeltaZpeak1Offset = FitResult.xindex(idx); 

        [FitResult.DeltaRpeak1, idx] = max(FitResult.RamideFit);
        FitResult.DeltaRpeak1Offset = FitResult.xindex(idx); 
        
    case 2
        % (2) fixed peak
        offs_amide = 3.5;
        [~, idx] = min(abs(Offset - offs_amide));
        FitResult.DeltaZpeak1 = FitResult.ZamideFit(idx);
        FitResult.DeltaRpeak1 = FitResult.RamideFit(idx);

        FitResult.DeltaZpeak1Offset = FitResult.xindex(idx); 
        FitResult.DeltaRpeak1Offset = FitResult.xindex(idx); 
end

% guan peak
FitParam.peak = 3;

[Ztmp, Rtmp] = CurveFunction(FitResult.Coefficents, FitResult.xindex, FitParam); 
FitResult.ZguanFit =FitResult.Background -Ztmp;
FitResult.RguanFit =Rtmp - Rbak;
switch flag
    case 1
        % (1) maximum peak
        [FitResult.DeltaZpeak2, idx] = max(FitResult.ZguanFit); 
        FitResult.DeltaZpeak2Offset = FitResult.xindex(idx); 

        [FitResult.DeltaRpeak2, idx] = max(FitResult.RguanFit);
        FitResult.DeltaRpeak2Offset = FitResult.xindex(idx); 
        
    case 2
        % (2) fixed peak
        offs_guan = 2;
        [~, idx] = min(abs(Offset - offs_guan));
        FitResult.DeltaZpeak2 = FitResult.ZguanFit(idx);
        FitResult.DeltaRpeak2 = FitResult.RguanFit(idx);

        FitResult.DeltaZpeak2Offset = FitResult.xindex(idx); 
        FitResult.DeltaRpeak2Offset = FitResult.xindex(idx); 
end

end