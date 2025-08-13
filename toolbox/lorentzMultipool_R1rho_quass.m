function y = lorentzMultipool_R1rho_quass(p, x, FitParam)
% Multi-pool Lorentzian function
% use R1rho signal equation: 
%       "Z = (1 - cos2thet*R1/R1rho).* exp(-R1rho * tsat) + cos2thet * R1/R1rho"
% Huabin Zhang   Email: huabinz@connect.hku.hk
% Input
%   p: exchange pool parameter, p(1) for Z0, p(2-4) for water, p(5-7) for pool_1, ...
%   x: frequency offset, [ppm]
%   FitParam: sequence-related parameters
%      FitParam.satpwr: saturation pulse power, [uT]
%      FitParam.Magfield: = B0 [T] * 42.5764 [Hz/T], [Hz]
%      FitParam.R1: = 1/T1, [Hz], 1 represents no T1 correction
%      FitParam.tsat: pulse duration, [s]
% Output
%   y: Z-spectra

satHz = FitParam.satpwr*42.5764;
cos2thet = 1 - satHz^2./(satHz^2+ (FitParam.Magfield*x).^2);

pn = (length(p)-1)/3;
R1rho = p(1);
dw_0 = p(4); % B0 shift / water peak shift
for i = 1:pn
    A = p(3*i-1);
    G = p(3*i);
    dw = p(3*i+1);
    if i == 1
        R1rho = R1rho + A*G^2/4./(G^2/4+(x-dw_0).^2);
    else
        R1rho = R1rho + A*G^2/4./(G^2/4+(x-dw_0-dw).^2);
    end
end

y = (1 - exp(- FitParam.R1 * FitParam.td) - cos2thet*FitParam.R1./R1rho).* exp(-R1rho * FitParam.tsat) + cos2thet * FitParam.R1./R1rho;

end