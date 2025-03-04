% 20241016
% simulate 3T Guan data, for PLOF quantitative comparison
% sequence parameters are consistent with HKU MR scanner
clear

addpath ..\toolbox\

%% load freq offs
offslist = "user1demo"; % m0 not included
offsFlag = 1;
fileID = fopen(offslist(offsFlag)+".txt", 'r');
data = cell2mat(textscan(fileID, '%f'));
fclose(fileID);
offs = reshape(data,[],1);

%% scanner settings
b0 = 3;
gamma = 267.5154109126009;
gamma_hz = gamma/2/pi;

%% saturation rf pulse 
pulse1_pwr = 0.8; % in uT
pulse1_dur = 2; % pulse duration in s
pulse_cell = {[pulse1_pwr*gamma_hz,0,pulse1_dur]};
pulse_tpost = 6.5e-3;

%% exchange settings
%       {name,            t1 [s],   t2 [s],    exch rate [Hz],  dw [ppm],    fraction (0~1)}
water  = {'water',        1.0,      0.04       1,               0,           1};
mt     = {'mt',           1.0,      4.0e-05,   30,              -2.5,        0.1351};
amide  = {'amide',        1.0,      0.1        50,              3.5,         300 * 0.0009009/100};
guanid = {'guanidine',    1.0,      0.1,       100,             2.0,         100 * 0.0009009/100};
noe    = {'noe',          1.3,      0.005,     20,              -3.5,        2000 * 0.0009009/100};

%% cest simulation
nf = length(offs);
zspec = zeros(nf,1); % all 5 pool
zspec_bak = zeros(nf,1); % water + mt + noe

tic
for idxoffs = 1:nf

    offs_temp = offs(idxoffs);
    % all 5 pools: water + mt + noe + amide + guan
    pools = {water; mt; amide; guanid; noe};
    magn = bmesolver(b0, gamma_hz, pools, pulse_cell, pulse_tpost, offs_temp, 0);
    zspec(idxoffs) = magn(length(pools)*2+1, end, end);

    % background: water + mt + noe
    pools = {water; mt; noe};
    magn = bmesolver(b0, gamma_hz, pools, pulse_cell, pulse_tpost, offs_temp, 0);
    zspec_bak(idxoffs) = magn(length(pools)*2+1, end, end);

end
figure();plot(offs,zspec,offs,zspec_bak)

%% save data
fileName = "zspec_"+offslist(offsFlag)+".mat";
save(fileName,"offs","zspec","zspec_bak",'-mat');
fprintf("data is saved to "+fileName+"\n")