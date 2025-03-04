function magn_end = bmesolver(magn_field, gamma_hz, pools, pulse_cell, pulse_tpost, offs, magn0)
% Bloch-Mcconnell Equations (BME) solver
% DATE: 2023/05/20
% AUTHOR: Jianpan Huang - jianpanhuang@outlook.com
%
% INPUT:
% magn_field: magnetic field in Tesla
% gamma: gyromagnetic ratio in MHz/T
% pools = {pool1; pool2; pool3; pool4;...}; pool1 is always free water pool
% pool = {name, T1, T2, exchange rate, chemical shift, proton fraction};
%         name: pool name
%         T1: T1 of the pool
%         T2: T2 of the pool
%         exchange rate: exchange rate of the exchangable proton, in Hz
%         chemicalShift: chemical shift of the exchangable proton, in ppm
%         proton fraction: proton fraction of the exchangable proton, 0～1
% pulseCell: {pulse1, pulse2, pulse3,...}
%     pulse: [power1,phase1,duration1;power2,phase2,duration2;...]
%     power: in Hz, 0 for no pulse but delay
%     phase: (0:x, pi/2:y, pi:-x, 3*pi/2:-y) 
%     duration: in s 
% pulse_repeat: number of repeat of the pulses in pulseCell
% pulse_tpost: post-pulse delay
% sequence: [pulse1, pulse2,...] -> [pulse1, pulse2,...] -> ... -> tpost -> acquisition
%           \_________________________________________________/
%                        repeat pulse_repeat times
% offs: frequency offset of the saturation pulse, in ppm
% magn0: initial magnetization vector in the follwoing order
%       [x(water); x(exch1); x(exch2);...
%       y(water); y(exch1); y(exch2);...
%       z(water); z(exch1); z(exch2);...]
%       z(water) = 1;
%       0: thermal equilibrium
%
% OUTPUT:
% magn_end: three dimensional matrix of magnetization vector after the pulse
%       dim 1: [x0; x1; x2; ...y0; y1; y2;...z0; z1; z2; ...]
%       dim 2: [after pulse 1; after pulse 2; after pulse 3;...]
%       dim 3: [after repeat 1; after repeat 2; after repeat 3;...]
%
% NOTE:
% The format of BME is the same as the one in JMR 154, 157–160 (2002) 

pools = pools(:);
water_pool = pools{1};
exch_pools = pools(2:end);
pool_num = length(pools);
pulse_num = length(pulse_cell);

%  fract: fractions for all pools
fract = cellfun(@(pool) pool{6},  pools);
%  magn0: initial magnetizations for all pools
%         [x0; x1; x2; ...y0; y1; y2;...z0; z1; z2; ...]
if (size(magn0(:)) ~= 3*pool_num)
    magn0 = [zeros(pool_num*2,1); fract(:)] / fract(1);
end

% create the relaxation matrix
r1 = cellfun(@(pool) 1/pool{2},  pools);
r2 = cellfun(@(pool) 1/pool{3},  pools);
r_mat = -diag([r2; r2; r1]);

% create the exchange matrix
exch2water = cellfun(@(pool)  pool{4},  exch_pools);
water2exch = cellfun(@(pool)  pool{6}*pool{4},  exch_pools) / water_pool{6};
exch_diag = [-sum(water2exch), -exch2water(:)'];
exch_block = diag(exch_diag);
exch_block(1,2:end) = exch2water(:)';
exch_block(2:end,1) = water2exch(:);
exch_mat = zeros(pool_num*3); % 3 dimensions
for ndims = 1:3
    exch_mat((ndims-1)*pool_num+1:ndims*pool_num, (ndims-1)*pool_num+1:ndims*pool_num) = exch_block;
end

% create the chemical shift (in ppm) matrix
dw_arr = cellfun(@(pool) pool{5},  pools); 
dw = 2*pi*magn_field*gamma_hz*(dw_arr-offs);        
cs_mat = diag([-dw; zeros(pool_num,1)],pool_num) ...
              + diag([dw; zeros(pool_num,1)], -pool_num); % chemical shift evolution frequency in radian

% create the evolution matrix for the pulse
for npulsenum = 1:pulse_num
    pulse = pulse_cell{npulsenum};
    for npulsestep = 1:size(pulse,1)
        xpul = ones(pool_num,1) * pulse(npulsestep,1) * cos(pulse(npulsestep,2)) * pi * 2;
        ypul = ones(pool_num,1) * pulse(npulsestep,1) * sin(pulse(npulsestep,2)) * pi * 2;
        zpul = zeros(pool_num,1);
        pulse_mat = diag([-zpul;-xpul], pool_num) + diag(ypul, pool_num*2) + ...
                    diag([zpul;xpul], -pool_num) + diag(-ypul, -pool_num*2);
        evol_mat = r_mat + exch_mat + cs_mat + pulse_mat;
        if npulsestep==1&&npulsenum==1
            magn_init = magn0;
        elseif npulsestep==1&&npulsenum~=1
            magn_init = magn_pulse(:,end,end);
        else
            magn_init = magn_pulsestep(:,npulsestep-1);
        end
        % Sove the equation using the formal solution: magn = (magn + ss_mat)*exp(evol_mat*t) - ss_mat
        % where ss_mat = evol_mat_inv*(-r_mat*magn0), times evol_mat on both sides
        % we get evol_mat*ss_mat = -r_mat*magn0, hence ss_mat = evol_mat/(-r_mat*magn0)
        t = pulse(npulsestep,3);
        ss_mat = mldivide(evol_mat, -r_mat*magn0);
        % ss_mat = inv(evol_mat)*(-r_mat)*magn0;
        magn_pulsestep(:,npulsestep) = expm(evol_mat*t)*(magn_init + ss_mat) - ss_mat;
    end
    magn_pulse(:,npulsenum) = magn_pulsestep(:,npulsestep);
end
magn_end(:,1) = magn_pulse(:,pulse_num);

% calculate the evolution after pulse if there is a tpost
if pulse_tpost ~= 0
    xpul = zeros(pool_num,1);
    ypul = zeros(pool_num,1);
    zpul = zeros(pool_num,1);
    pulse_mat = diag([-zpul;-xpul], pool_num) + diag(ypul, pool_num*2) + ...
                diag([zpul;xpul], -pool_num) + diag(-ypul, -pool_num*2);
    evol_mat = r_mat + exch_mat + cs_mat + pulse_mat;
    magn_init = magn_end;
    % Sove the equation using the formal solution: magn = (magn + ss_mat)*exp(evol_mat*t) - ss_mat
    % where ss_mat = evol_mat_inv*(-r_mat*magn0), times evol_mat on both sides
    % we get evol_mat*ss_mat = -r_mat*magn0, hence ss_mat = evol_mat/(-r_mat*magn0)
    t = pulse_tpost;
    ss_mat = mldivide(evol_mat, -r_mat*magn0);
    % ss_mat = inv(evol_mat)*(-r_mat)*magn0;
    magn_end(:,1) = expm(evol_mat*t)*(magn_init + ss_mat) - ss_mat;
end
end