% 20241016
% simulate 3T Guan data, for PLOF quantitative comparison
% sequence parameters are consistent with HKU MR scanner
clear

% ----- Cross-platform path handling -----
here = fileparts(mfilename('fullpath'));  % folder of this script
if isempty(here)                          % e.g., if run as a Live Script
    here = pwd;
end

% Add toolbox folder (.. / toolbox)
toolboxDir = fullfile(here, '..', 'toolbox');
if exist(toolboxDir, 'dir')
    addpath(toolboxDir);
else
    warning('Toolbox folder not found: %s', toolboxDir);
end

%% load freq offs
offslist = "user1demo"; % m0 not included
offsFlag = 1;

% If the .txt is next to this script:
offsFile = fullfile(here, offslist(offsFlag) + ".txt");

% If it lives somewhere else, e.g., in a "data" subfolder:
% offsFile = fullfile(here, 'data', offslist(offsFlag) + ".txt");

fileID = fopen(char(offsFile), 'r');
if fileID < 0
    error('Failed to open offsets file: %s', offsFile);
end
data = cell2mat(textscan(fileID, '%f'));
fclose(fileID);
offs = reshape(data, [], 1);

%% scanner settings
b0 = 3;
gamma = 267.5154109126009;
gamma_hz = gamma/2/pi;

%% saturation rf pulse
pulse1_pwr = 0.8; % in uT
pulse1_dur = 0.5; % pulse duration in s
pulse_cell = {[pulse1_pwr*gamma_hz, 0, pulse1_dur]};
pulse_tpost = 6.5e-3;
t_delay = pulse1_dur; % relaxation delay

%% exchange settings
%       {name,            t1 [s],   t2 [s],    exch rate [Hz],  dw [ppm],    fraction (0~1)}
water  = {'water',        1.0,      0.04       1,               0,           1};
mt     = {'mt',           1.0,      4.0e-05,   30,              -2.5,        0.1351};
amide  = {'amide',        1.0,      0.1        50,              3.5,         300 * 0.0009009/100};
guanid = {'guanidine',    1.0,      0.1,       100,             2.0,         100 * 0.0009009/100};
noe    = {'noe',          1.3,      0.005,     20,              -3.5,        2000 * 0.0009009/100};

%% cest simulation
tic
offset_ref = 300; % ppm, reference scan
% background: water + mt + noe + amide + guan
pools = {water; mt; amide; guanid; noe};
zspec = bmesimulator(offs, t_delay, pools, pulse_cell, pulse_tpost, offset_ref, b0, gamma_hz);
% background: water + mt + noe
pools = {water; mt; noe};
zspec_bak = bmesimulator(offs, t_delay, pools, pulse_cell, pulse_tpost, offset_ref, b0, gamma_hz);
figure(); plot(offs, zspec, offs, zspec_bak)

%% save data (to the script folder or a subfolder)
outDir = here;  % or: fullfile(here, 'output')
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
fileName = "zspec_" + offslist(offsFlag) + ".mat";
outPath = fullfile(outDir, fileName);
save(char(outPath), "offs", "zspec", "zspec_bak", '-mat');
fprintf("data is saved to %s\n", outPath);


function zspec = bmesimulator(offs, t_delay, pools, pulse_cell, pulse_tpost, offset_ref, B0, gamma_hz)
    
    % ---------- reference scan Z -------------
    poolnum = size(pools);
    
    % % same effect as just a exponential model
    lastItems = cellfun(@(pool) pool{end}, pools);  % e.g., [1, 0.001]
    magn_zeros = zeros([poolnum(1)*3 1]);
    magn_zeros(poolnum(1)*2+1:end) = lastItems;
    magn_rec = magn_zeros;
    disp(magn_rec)
    for ref_ss=1:3
        magn = bmesolver(B0, gamma_hz, pools, pulse_cell, pulse_tpost, offset_ref, magn_rec);
        fprintf('%d th TR for steady state:\n',ref_ss)
        disp(magn)
        magn_spoil = magn*0; 
        magn_spoil(1:length(pools)) = magn(end-length(pools)+1:end); 
        z_ref = magn(length(pools)*2+1, end, end);
        disp(magn_spoil)
        disp(z_ref)
        magn_rec = bmesolver(B0, gamma_hz, pools, [], t_delay, 0, magn_spoil);
    end
    
    % ----------- simulate Z spectrum --------------
    nf = length(offs);
    zspec = zeros(nf,1);
    for idxoffs = 1:nf
        offs_temp = offs(idxoffs);
    
        % (1) CEST saturation, R1rho decay
        magn = bmesolver(B0, gamma_hz, pools, pulse_cell, pulse_tpost, offs_temp, magn_rec);
    
        % (2) FA = 90, Mz --> Mx
        magn_spoil = magn*0; 
        magn_spoil(1:length(pools)) = magn(end-length(pools)+1:end); % actually no impact because Mxy is spoiled after saturation
    
        % (3) readout, single proportional to Mxy
        zspec(idxoffs) = magn(length(pools)*2+1, end, end);
    
        % (4) R1 relaxation
        magn_rec = bmesolver(B0, gamma_hz, pools, [], t_delay, 0, magn_spoil);
    end

    zspec = zspec/z_ref;
end