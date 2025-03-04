function [Z, R] = CurveFunction(x, xdata, FitParam)
%------------ Input------------%
%
%   x is the list of Lorentzian-peak related fitting parameters: 
%       x(1~4) for background;   
%       x(5~7) for peak_1;   x(8~10) for peak_2; ...
%
%   x =  [Amp_background, width_background, Coefficient_of_Zero_Order_background, Coefficient_of_First_Order_background, ...
%           Rexch_peak_1, width_peak_1, center_peak_1, ...
%           Rexch_peak_2, width_peak_2, center_peak_2, ...]
%
%   xdata is the list of offsets, y value (signal Z) is in the FitParam as
%   FitParam.Satruration
%
%   peak = FitParam.peak;
%   peak = 0: only fit background;
%           = 1: all peaks are used for the fitting;
%           = 2: only backgound and peak_1 are used for the fitting;
%           = 3: only backgound and peak_2 are used for the fitting;
%           ...
%
%------------ Output------------%
% 
%   Z is the Z-spectrum.


% Eq. (6): R_back = ; C_0/1/2/3 becomes x(4/5/6/7)
peakoffset = FitParam.PeakOffset;
switch length(x)
    case 10 % first-order polynomial
        % x_poly = x(1:4);
        Rbak = x(1)*x(2)^2./( x(2)^2 + 4*(xdata ).^2 ) + x(3) + x(4)*0.001*(xdata - peakoffset);
        x_cest = x(5:10);
    case 11 % second-order polynomial
        % x_poly = x(1:5);
        Rbak = x(1)*x(2)^2./( x(2)^2 + 4*(xdata ).^2 ) + x(3) + x(4)*0.001*(xdata - peakoffset)...
                + x(5)*0.001*(xdata - peakoffset).^2;
        x_cest = x(6:11);
end

 
% Eq. (5): R = R_exch * (w/2) ^2 over (w/2) ^2 + (delta_w - delta_w_exch)^2
Ramide = x_cest(1)*x_cest(2)^2./( x_cest(2)^2 + 4*(xdata - x_cest(3)).^2 );
Rguan = x_cest(4)*x_cest(5)^2./( x_cest(5)^2 + 4*(xdata - x_cest(6)).^2 );

 % cos^2(theta)
satHz = FitParam.satpwr*42.58;
cos2thet = 1 - satHz^2./(satHz^2+ (FitParam.Magfield*xdata).^2);

peak = FitParam.peak;
if peak == 0 
    Rall = Rbak;
elseif peak == 1 
    Rall = Rbak + Ramide + Rguan;
elseif peak == 2
    Rall = Rbak + Ramide;
elseif peak == 3
    Rall = Rbak + Rguan;
elseif peak == 21
    Rall = x(3) + Ramide;
elseif peak == 31
    Rall = x(3) + Rguan;
end

% Output of Z-spectrum 
Z = (1 - cos2thet*FitParam.R1./Rall).* exp(-Rall * FitParam.tsat) + cos2thet * FitParam.R1./Rall;
R = Rall;
end