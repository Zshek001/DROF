% 20250121
% simulate 3T Guan data, for PLOF quantitative comparison
% sequence parameters are consistent with HKU MR scanner
%   using random data set
clear
Npixel = 5000; % Zspec number
rng default % for reproducibility

%% load freq offs
fileID = fopen("user19SSFSE.txt", 'r');
data = cell2mat(textscan(fileID, '%f'));
fclose(fileID);
offs = reshape(data,[],1);
nf = length(offs);

%% scanner settings
b0 = 3;
gamma = 267.5154109126009;
gamma_hz = gamma/2/pi;

%% saturation rf pulse 
pulse1_pwr = 0.8; % in uT
pulse1_dur = 2; % pulse duration in s
pulse_cell = {[pulse1_pwr*gamma_hz,0,pulse1_dur]};
pulse_tpost = 6.5e-3;

%% simulation
paraList = zeros(7,Npixel);
paraInfo = ["MTf_relative","Guank_Hz","Amidef_mM","Guanf_mM","NOEf_mM","Amidek_Hz","NOEk_Hz"];
zspec = zeros(nf,Npixel); % all 5 pool
zspec_bak = zeros(nf,Npixel); % water + mt + noe

tic
backNum = 0;
for idxPixel = 1:Npixel
    fprintf(repmat('\b',1,backNum));
    backNum =fprintf('Calculation process: %d/%d',idxPixel,Npixel);
    % mt
    a = 5000; b = 25000; fmt_mM = (b-a).*rand(1,1) + a; % mM

    % amide
    a = 20; b = 400; famide_mM = (b-a).*rand(1,1) + a; % mM

    % guan
    a = 10; b = 100; fguan_mM = (b-a).*rand(1,1) + a; % mM
    a = 50; b = 400; kguan = (b-a).*rand(1,1) + a;

    % noe
    a = 100; b = 2000; fnoe_mM = (b-a).*rand(1,1) + a; % mM
    
    %% exchange settings, water proton 111M
    %       {name,            t1 [s],   t2 [s],    exch rate [Hz],  dw [ppm],    fraction (0~1)}
    water  = {'water',        1,        0.04,       1,               0,           1};
    mt     = {'mt',           1.0,      4.0e-05,    30,             -2.5,         fmt_mM * 0.0009009/100};
    amide  = {'amide',        1.0,      0.1,        50,              3.5,         famide_mM * 0.0009009/100};
    guanid = {'guanidine',    1.0,      0.1,        kguan,           2.0,         fguan_mM * 0.0009009/100};
    noe    = {'noe',          1.3,      0.001,      25,             -3.5,         fnoe_mM * 0.0009009/100};

    paraList(1,idxPixel) = fmt_mM;
    paraList(2,idxPixel) = kguan;
    paraList(3,idxPixel) = famide_mM;
    paraList(4,idxPixel) = fguan_mM;
    paraList(5,idxPixel) = fnoe_mM;


    %% cest simulation
    % Z-spectrum
    for idxoffs = 1:nf
        offs_temp = offs(idxoffs);

        % all 5 pools: water + mt + noe + amide + guan
        pools = {water; mt; amide; guanid; noe};
        magn = bmesolver(b0, gamma_hz, pools, pulse_cell, pulse_tpost, offs_temp, 0);
        zspec(idxoffs,idxPixel) = magn(length(pools)*2+1, end, end);

        % background: water + mt + noe
        pools = {water; mt; noe};
        magn = bmesolver(b0, gamma_hz, pools, pulse_cell, pulse_tpost, offs_temp, 0);
        zspec_bak(idxoffs,idxPixel) = magn(length(pools)*2+1, end, end);
    end
end 
fprintf("\nDone\n")
fprintf('elapsed time: %.4f s\n', toc);

%% save data
fileName = "simData_random"+num2str(Npixel)+".mat";
save(fileName,"offs","zspec","zspec_bak","paraList","paraInfo",'-mat');
fprintf("data is saved to "+fileName+"\n")