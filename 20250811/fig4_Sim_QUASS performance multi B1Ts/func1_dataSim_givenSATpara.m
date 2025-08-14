function [zspecList,paraList] = func1_dataSim_givenSATpara(SATpara, offs, Npixel)
% simulate random zspec at given SATpara
% INPUT
%   SATpara: [1,3], [B1(uT), tsat(s), tdelay(s)]
%   offs: [1,nf], including M0, three 200 ppm and one -300 ppm
% OUTPUT
%   zspecList: [nf,Npixel]
%   paraList: [4,Npixel], concentration of 4 pools

%% simulation with TR
nf = length(offs);
zspecList = zeros(nf,Npixel);
paraList = zeros(4,Npixel);

pulse1_pwr = SATpara(1); % in uT
pulse1_dur = SATpara(2); % pulse duration in s, ts
Td = SATpara(3); % tdelay or trecovery in s, td

for idxPixel = 1:Npixel
    % (1) scanner settings
    B0 = 3;
    gamma = 267.5154109126009;
    gamma_hz = gamma/2/pi;
    
    % (2) saturation rf pulse 
    pulse1_phase = zeros(size(pulse1_dur)); % pulse duration in s
    pulse1 = [pulse1_pwr*gamma_hz, pulse1_phase, pulse1_dur];
    pulse2_tdelay = 0.00; % pulse delay in s
    pulse2 = [0, 0, pulse2_tdelay];
    pulse_block = {pulse1, pulse2};
    pulse_block_len = length(pulse_block);
    pulse_repeat = 1; % repeat number of pulse cell
    pulse_cell=cell(1,pulse_block_len*pulse_repeat);
    for idxoffs=1:pulse_block_len:size(pulse_cell,2)
        pulse_cell(1,idxoffs:idxoffs+pulse_block_len-1) = pulse_block;
    end
    last_delay = 0; % 1 ends with delay, 0 ends with pulse (no last delay)
    if last_delay == 0
        pulse_cell(end) = [];
    end
    pulse_tpost = 0;
    
    % (3) exchange settings
    % mt
    a = 5000; b = 25000; fmt_mM = (b-a).*rand(1,1) + a; % mM
    a = 20; b = 60; kMT = (b-a).*rand(1,1) + a;

    % amide
    a = 20; b = 400; famide_mM = (b-a).*rand(1,1) + a; % mM
    a = 20; b = 50; kamide = (b-a).*rand(1,1) + a;

    % guan
    a = 10; b = 100; fguan_mM = (b-a).*rand(1,1) + a; % mM
    a = 200; b = 500; kguan = (b-a).*rand(1,1) + a;

    % noe
    a = 100; b = 2000; fnoe_mM = (b-a).*rand(1,1) + a; % mM
    a = 10; b = 30; knoe = (b-a).*rand(1,1) + a;

    %       {name,            t1 [s],   t2 [s],    exch rate [Hz],  dw [ppm],    fraction (0~1)}
    water  = {'water',        1,        0.04,       1,               0,           1};
    mt     = {'mt',           1.0,      4.0e-05,    kMT,             -2.5,         fmt_mM * 0.0009009/100};
    amide  = {'amide',        1.0,      0.1,        kamide,              3.5,         famide_mM * 0.0009009/100};
    guanid = {'guanidine',    1.0,      0.1,        kguan,           2.0,         fguan_mM * 0.0009009/100};
    noe    = {'noe',          1.3,      0.001,      knoe,             -3.5,         fnoe_mM * 0.0009009/100};
    pools = {water; mt; amide; guanid; noe};

    paraList(:,idxPixel) = [fmt_mM;famide_mM;fguan_mM;fnoe_mM];
    
    % (4) BME solver for Z-spectrum
    zspec = zeros(nf,1);
    for n = 1:nf
        offs_temp = offs(n);

        % 1) CEST saturation, R1rho decay
        if n == 1
            magn = bmesolver(B0, gamma_hz, pools, pulse_cell, pulse_tpost, offs_temp, 0);
        else
            magn = bmesolver(B0, gamma_hz, pools, pulse_cell, pulse_tpost, offs_temp, magn_rec);
        end
        % 2) FA = 90, Mz --> Mx
        magn_spoil = magn*0; 
        magn_spoil(1:length(pools)) = magn(end-length(pools)+1:end); % actually no impact because Mxy is spoiled after saturation

        % 3) readout, signale propotional to Mxy
        zspec(n,1) = magn(length(pools)*2+1, end, end);

        % 4) R1 relaxation
        magn_rec = bmesolver(B0, gamma_hz, pools, [], Td, 0, magn_spoil);
        % magn_rec(end:end-3) = 0;
        % magn_rec = 1 - exp(-1*Td);
        
    end

    zspecList(:,idxPixel) = zspec;
end

%% truncate M0
zspecM0 = zspecList(4,:);
zspecList = zspecList./zspecM0;
zspecList(1:4,:) = [];
offs(1:4) = [];

end