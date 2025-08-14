function zspec = zspecModel_QUASS_m0us(par, offs_Hz)
% Zspec in QUASS format
% Jianpan Huang   Email: jianpanhuang@outlook.com
% Ref: Sun Philip-2021. MRM 86(2): 765-776.
%   for control scan (e.g., -300 ppm), assum unsteady-state, i.e. Z(-300ppm)=1

R1rho = par(1); % Hz
B1 = par(2); % uT
Ts = par(3); % s
Td = par(4); % s
R1w = par(5); % Hz

gamma = 267.5154109126009/2/pi; % for protons, in 1e6 Hz/T = 1 Hz/uT
theta = atan(B1*gamma./offs_Hz); % tilted angle for effective field

zspec = (1-exp(-R1w*Td)) * exp(-R1rho*Ts) + ...
    R1w/R1rho*cos(theta).^2*(1-exp(-R1rho*Ts));
zspec = zspec ./ (1-exp(-R1w*(Ts+Td)));