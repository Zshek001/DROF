% 20241218
% simulate zspec in CSF/GM/WM based on the pre-calculated mask
%   with gaussian noise,
clear
% close all

%% load brain mask
load('roi_csf_gm_wm.mat','roi_csf','roi_gm','roi_wm');
% Fig1 = figure();set(gcf,'Position',[150 150 800 250]);
% tiledlayout(1,3,"TileSpacing","compact","Padding","compact")
% nexttile;imshow(roi_csf,[]);title("ROI csf");
% nexttile;imshow(roi_gm,[]);title("ROI gm");
% nexttile;imshow(roi_wm,[]);title("ROI wm");

%% load freq offs
load('roi_csf_gm_wm.mat','offs'); % include -300
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

%% exchange pools of CSF, GM, WM
% WM
%       {name,            t1 [s],   t2 [s],    exch rate [Hz],  dw [ppm],    fraction (0~1)}
water  = {'water',        1.05,     0.099,      1,               0,           0.99};
mt     = {'mt',           1.05,     4.4e-05,    29,             -3,           15500 * 0.0009009/100};
amide  = {'amide',        1.3,      0.0013,     40,              3.5,         350 * 0.0009009/100};
guanid = {'guanidine',    1.3,      0.009,      400,             2,           75 * 0.0009009/100};
PCr    = {'PCr',          1.6,      0.01,       300,             2.5,         5 * 0.0009009/100};
noe    = {'noe',          1.3,      0.00045,     55,             -3.5,         1800 * 0.0009009/100};
noe2   = {'noe2',         1.6,      0.00045,     20,             -1.8,         100 * 0.0009009/100};
amine  = {'amine',        1.3,      0.01,       2000,            3,           40 * 0.0009009/100};
pools_wm = {water; mt; amide; guanid; noe;amine;noe2;PCr};

% GM
%       {name,            t1 [s],   t2 [s],    exch rate [Hz],  dw [ppm],    fraction (0~1)}
water  = {'water',        1.9,      0.25,       1,               0,           0.99};
mt     = {'mt',           1.6,      5.3e-05,    39,             -2.5,         5000 * 0.0009009/100};
amide  = {'amide',        1.6,      0.005,      60,              3.5,         110 * 0.0009009/100};
guanid = {'guanidine',    1.6,      0.03,       400,             1.9,         50 * 0.0009009/100};
PCr    = {'PCr',          1.6,      0.035,      150,             2.5,         45 * 0.0009009/100};
noe    = {'noe',          1.6,      0.0035,     55,             -3.5,         180 * 0.0009009/100};
noe2    = {'noe2',        1.6,      0.003,      20,             -1.8,         400 * 0.0009009/100};
amine  = {'amine',        1.6,      0.005,      1000,            2.9,         30 * 0.0009009/100};
pools_gm = {water; mt; amide; guanid; noe;amine;noe2;PCr};

% CSF
%       {name,            t1 [s],   t2 [s],    exch rate [Hz],  dw [ppm],    fraction (0~1)}
water  = {'water',        3.5,      1.5,         1,             0,           0.99};
mt     = {'mt',           3.5,      3.7e-05,     29,           -3,           200 * 0.0009009/100};
amide  = {'amide',        3.5,      0.001,       50,            3.5,         40 * 0.0009009/100};
guanid = {'guanidine',    3.5,      0.009,       400,           2,           0 * 0.0009009/100};
PCr    = {'PCr',          1.6,      0.035,       300,           2.5,         0 * 0.0009009/100};
noe    = {'noe',          3.5,      0.0005,      16,           -3.5,         100 * 0.0009009/100};
noe2    = {'noe2',        1.6,      0.003,       20,           -1.8,         0 * 0.0009009/100};
amine  = {'amine',        3.5,      0.01,        3000,          3,           5 * 0.0009009/100};
pools_csf = {water; mt; amide; guanid; noe;amine;noe2;PCr};

%% simulation
zspec = zeros(nf,3); % 5 pools, 3 roi [csf, gm, wm]
for idxoffs = 1:nf
    offs_temp = offs(idxoffs);

    % all 5 pools: water + mt + noe + amide + guan
    magn = bmesolver(b0, gamma_hz, pools_csf, pulse_cell, pulse_tpost, offs_temp, 0);
    zspec(idxoffs,1) = magn(length(pools_csf)*2+1, end, end);

    magn = bmesolver(b0, gamma_hz, pools_gm, pulse_cell, pulse_tpost, offs_temp, 0);
    zspec(idxoffs,2) = magn(length(pools_gm)*2+1, end, end);

    magn = bmesolver(b0, gamma_hz, pools_wm, pulse_cell, pulse_tpost, offs_temp, 0);
    zspec(idxoffs,3) = magn(length(pools_wm)*2+1, end, end);
end

% m0 normalize and truncate
zspec = zspec./zspec(1,:); % normalize using m0 = -300ppm
zspec = zspec(2:end,:); % exclude -300
offs = offs(2:end);

%% construct digital brain
rng default % for reproducibility
sigma = 0.002;

zspec_brain = zeros([size(roi_csf),length(offs)]);
for i = 1:length(offs)
    img_2d = zeros(size(roi_csf));
    img_2d(roi_csf==1) = zspec(i,1) + sigma * randn(1,length(find(roi_csf==1)));
    img_2d(roi_gm==1) = zspec(i,2) + sigma * randn(1,length(find(roi_gm==1)));
    img_2d(roi_wm==1) = zspec(i,3) + sigma * randn(1,length(find(roi_wm==1)));

    zspec_brain(:,:,i) = img_2d;
end

%% 1D zspec comparision
offssim = offs;
pixelLoc_WM = [104,152];
pixelLoc_GM = [215,178];
pixelLoc_CSF = [156,152];
load('.\20250203-invivoresult\output_DZ_LDA.mat','offs','zmap_mea');

Fig0 = figure();set(gcf,'Position',[50 50 400 300]);
hold on
plot(offs,squeeze(zmap_mea(pixelLoc_WM(2),pixelLoc_WM(1),:)),'r',offssim,squeeze(zspec_brain(pixelLoc_WM(2),pixelLoc_WM(1),:)),'r--'); %WM
plot(offs,squeeze(zmap_mea(pixelLoc_GM(2),pixelLoc_GM(1),:)),'b',offssim,squeeze(zspec_brain(pixelLoc_GM(2),pixelLoc_GM(1),:)),'b--'); %GM
plot(offs,squeeze(zmap_mea(pixelLoc_CSF(2),pixelLoc_CSF(1),:)),'k',offssim,squeeze(zspec_brain(pixelLoc_CSF(2),pixelLoc_CSF(1),:)),'k--'); %CSF
hold off;set(gca,'XDir','reverse');
xlim([-6,6])
legend('WM (mea)','WM (sim)','GM (mea)','GM (sim)','CSF (mea)','CSF (sim)','Location','southeast');
title('simulation vs. measurement')

%% save
save('digital_brain_wnoise.mat',"zspec_brain",'offs','roi_csf','roi_gm','roi_wm');

%% CEST analysis
clear

addpath ..\toolbox\
addpath ..\toolbox\fitting_main\
multistart_N = 5; % times for MultiStart

OutputDir = ".\Outprocess";
DirTemp = OutputDir;
cnt = 1;
while(1)
    if exist(DirTemp,'dir')
        DirTemp = OutputDir + "_" + num2str(cnt);
        cnt = cnt + 1;
    else
        break;
    end
end
OutputDir = DirTemp + "\";
mkdir(OutputDir);

load('digital_brain_wnoise.mat',"zspec_brain",'offs','roi_csf','roi_gm','roi_wm');
zmap_mea = squeeze(abs(zspec_brain)); % [nx,ny,nf]
[nx,ny,nf] = size(zmap_mea);

roi = logical(roi_csf+roi_gm+roi_wm);
indnz_csf = find(roi_csf~=0);
indnz_gm = find(roi_gm~=0);
indnz_wm = find(roi_wm~=0);

% indnz = [indnz_csf(1),indnz_gm(1),indnz_wm(1)]; % 3 pixels
indnz = [indnz_csf;indnz_gm;indnz_wm]; % roi pixels

indnz_vec_csf = 1:length(indnz_csf);
indnz_vec_gm = 1+length(indnz_csf):length(indnz_csf)+length(indnz_gm);
indnz_vec_wm = 1+length(indnz_csf)+length(indnz_gm):length(indnz);

zspec_vec = reshape(permute(zmap_mea,[3,1,2]),nf,[]); % [nf,Npixel]
zspec_nz_vec = zspec_vec(:,indnz);

% to be output
amide_DZ_map = zeros(nx,ny);
guan_DZ_map = zeros(nx,ny);
noe_DZ_map = zeros(nx,ny);
amide_CO_map = zeros(nx,ny);
guan_CO_map = zeros(nx,ny);
noe_CO_map = zeros(nx,ny);

% (1) Delta Z
processing_DZ_methods = {@fitting_MTR,@fitting_LDA};
output_filenames = {'DZ_MTR','DZ_LDA'};
for i = 1:length(processing_DZ_methods)
    method = processing_DZ_methods{i};
    [Zamide_DZ_vec, Zguan_DZ_vec, Znoe_DZ_vec] = method(offs, zspec_nz_vec, multistart_N);

    amide_DZ_map(indnz_csf) = Zamide_DZ_vec(indnz_vec_csf); amide_DZ_map(indnz_gm) = Zamide_DZ_vec(indnz_vec_gm); amide_DZ_map(indnz_wm) = Zamide_DZ_vec(indnz_vec_wm);
    guan_DZ_map(indnz_csf) = Zguan_DZ_vec(indnz_vec_csf); guan_DZ_map(indnz_gm) = Zguan_DZ_vec(indnz_vec_gm); guan_DZ_map(indnz_wm) = Zguan_DZ_vec(indnz_vec_wm);
    noe_DZ_map(indnz_csf) = Znoe_DZ_vec(indnz_vec_csf); noe_DZ_map(indnz_gm) = Znoe_DZ_vec(indnz_vec_gm); noe_DZ_map(indnz_wm) = Znoe_DZ_vec(indnz_vec_wm);

    save(OutputDir+"output_"+output_filenames{i}+".mat",'amide_DZ_map','guan_DZ_map','noe_DZ_map','roi','offs','zmap_mea');

    fprintf('saved to %s\n', "output_"+output_filenames{i}+".mat");
end

% (2) Coefficience, MPLF
processing_DZ_methods = {@fitting_MPLF};
output_filenames = {'CO_MPLF'};
for i = 1:1
    method = processing_DZ_methods{1};
    [Zamide_DZ_vec, Zguan_DZ_vec, Znoe_DZ_vec,fit_paratemp] = method(offs, zspec_nz_vec, multistart_N);

    amide_DZ_map(indnz_csf) = Zamide_DZ_vec(indnz_vec_csf); amide_DZ_map(indnz_gm) = Zamide_DZ_vec(indnz_vec_gm); amide_DZ_map(indnz_wm) = Zamide_DZ_vec(indnz_vec_wm);
    guan_DZ_map(indnz_csf) = Zguan_DZ_vec(indnz_vec_csf); guan_DZ_map(indnz_gm) = Zguan_DZ_vec(indnz_vec_gm); guan_DZ_map(indnz_wm) = Zguan_DZ_vec(indnz_vec_wm);
    noe_DZ_map(indnz_csf) = Znoe_DZ_vec(indnz_vec_csf); noe_DZ_map(indnz_gm) = Znoe_DZ_vec(indnz_vec_gm); noe_DZ_map(indnz_wm) = Znoe_DZ_vec(indnz_vec_wm);

    amide_CO_map(indnz_csf) = fit_paratemp(5,indnz_vec_csf); amide_CO_map(indnz_gm) = fit_paratemp(5,indnz_vec_gm); amide_CO_map(indnz_wm) = fit_paratemp(5,indnz_vec_wm);
    guan_CO_map(indnz_csf) = fit_paratemp(14,indnz_vec_csf); guan_CO_map(indnz_gm) = fit_paratemp(14,indnz_vec_gm); guan_CO_map(indnz_wm) = fit_paratemp(14,indnz_vec_wm);
    noe_CO_map(indnz_csf) = fit_paratemp(8,indnz_vec_csf); noe_CO_map(indnz_gm) = fit_paratemp(8,indnz_vec_gm); noe_CO_map(indnz_wm) = fit_paratemp(8,indnz_vec_wm);

    nterm = size(fit_paratemp,1);
    fit_para = zeros(nterm,nx*ny);
    fit_para(:,indnz_csf) = fit_paratemp(:,indnz_vec_csf);
    fit_para(:,indnz_gm) = fit_paratemp(:,indnz_vec_gm);
    fit_para(:,indnz_wm) = fit_paratemp(:,indnz_vec_wm);
    fit_para = reshape(fit_para,[nterm,nx,ny]);

    save(OutputDir+"output_"+output_filenames{1}+".mat",'amide_DZ_map','guan_DZ_map','noe_DZ_map','roi','offs','fit_para',...
                                                                    'amide_CO_map','guan_CO_map','noe_CO_map');
    fprintf('saved to %s\n', "output_"+output_filenames{1}+".mat");
end

% (3) DROF
B0 = 3; % T
FitParam.R1 = 1; % brain ~ 1s @ 3T
FitParam.satpwr = 0.8; % [uT]
FitParam.tsat = 2; % [s]
FitParam.Magfield = 42.5764 * B0; % [Hz]
FitParam.PeakOffset = 3.5; % used in PLOF, position to be polynomialize

processing_R1rho_methods = {@fitting_DROF};
output_filenames = {'DROF'};
for i = 1:length(processing_R1rho_methods)
    method = processing_R1rho_methods{i};
    [Zamide_DZ_vec, Zguan_DZ_vec, Znoe_DZ_vec, fit_paratemp] = method(offs, zspec_nz_vec, FitParam, multistart_N); 

    amide_DZ_map(indnz_csf) = Zamide_DZ_vec(indnz_vec_csf); amide_DZ_map(indnz_gm) = Zamide_DZ_vec(indnz_vec_gm); amide_DZ_map(indnz_wm) = Zamide_DZ_vec(indnz_vec_wm);
    guan_DZ_map(indnz_csf) = Zguan_DZ_vec(indnz_vec_csf); guan_DZ_map(indnz_gm) = Zguan_DZ_vec(indnz_vec_gm); guan_DZ_map(indnz_wm) = Zguan_DZ_vec(indnz_vec_wm);
    noe_DZ_map(indnz_csf) = Znoe_DZ_vec(indnz_vec_csf); noe_DZ_map(indnz_gm) = Znoe_DZ_vec(indnz_vec_gm); noe_DZ_map(indnz_wm) = Znoe_DZ_vec(indnz_vec_wm);

    amide_CO_map(indnz_csf) = fit_paratemp(5,indnz_vec_csf); amide_CO_map(indnz_gm) = fit_paratemp(5,indnz_vec_gm); amide_CO_map(indnz_wm) = fit_paratemp(5,indnz_vec_wm);
    guan_CO_map(indnz_csf) = fit_paratemp(14,indnz_vec_csf); guan_CO_map(indnz_gm) = fit_paratemp(14,indnz_vec_gm); guan_CO_map(indnz_wm) = fit_paratemp(14,indnz_vec_wm);
    noe_CO_map(indnz_csf) = fit_paratemp(8,indnz_vec_csf); noe_CO_map(indnz_gm) = fit_paratemp(8,indnz_vec_gm); noe_CO_map(indnz_wm) = fit_paratemp(8,indnz_vec_wm);

    nterm = size(fit_paratemp,1);
    fit_para = zeros(nterm,nx*ny);
    fit_para(:,indnz_csf) = fit_paratemp(:,indnz_vec_csf);
    fit_para(:,indnz_gm) = fit_paratemp(:,indnz_vec_gm);
    fit_para(:,indnz_wm) = fit_paratemp(:,indnz_vec_wm);
    fit_para = reshape(fit_para,[nterm,nx,ny]);

    save(OutputDir+"output_"+output_filenames{1}+".mat",'amide_DZ_map','guan_DZ_map','noe_DZ_map','roi','offs','fit_para',...
                                                                    'amide_CO_map','guan_CO_map','noe_CO_map');
    fprintf('saved to %s\n', "output_"+output_filenames{i}+".mat");
end

%% display
load('20241029-colormap-black-bg.mat','colormap_Jet_black_bg','colormap_Parula_black_bg','colormap_Turbo_black_bg','colormap_Inferno_bbg')
% load data
filepath = ".";
filenameList = ["output_DZ_MTR.mat","output_DZ_LDA.mat","output_CO_MPLF.mat","output_DROF.mat"];
% store variables:
%   pools_CO_map
%   pools_DZ_map
%   fit_para
%   offs, roi, zmap_mea, img_m0

amide_map = [];
guan_map = [];
noe_map = [];
mt_map = [];
for i = 1:4
    load(OutputDir+filepath+"\"+filenameList(i));
    if i <= 2
        if i == 1 % MTR background
            amide_DZ_map(roi==0) = -100;
            guan_DZ_map(roi==0) = -100;
            noe_DZ_map(roi==0) = -100;
        end
        amide_map = cat(3,amide_map,amide_DZ_map);
        guan_map = cat(3,guan_map,guan_DZ_map);
        noe_map = cat(3,noe_map,noe_DZ_map);
    elseif i == 3 % MPLF, Z [%]
        amide_map = cat(3,amide_map,amide_CO_map*100);
        guan_map = cat(3,guan_map,guan_CO_map*100);
        noe_map = cat(3,noe_map,noe_CO_map*100);
        mt_map = cat(3,mt_map,squeeze(fit_para(11,:,:))*100);
    elseif i == 4 % DROF, R1rho [ms^{-1}]
        amide_map = cat(3,amide_map,amide_CO_map*1000);
        guan_map = cat(3,guan_map,guan_CO_map*1000);
        noe_map = cat(3,noe_map,noe_CO_map*1000);
        mt_map = cat(3,mt_map,squeeze(fit_para(11,:,:))*1000);
    end
end

% display
tag = ["Z [%]","Z [%]","Z [%]","R1rho [ms^{-1}]"];
amide_ran = [[-6,2]; [0,25]; [0,7]; [0,150]];
guan_ran = [[-2,0.5]; [0,16]; [0,10]; [0,75]];
noe_ran = [[-1,6]; [0,23]; [0,13.5]; [0,230]];
mt_ran = [[0,21]; [0,310]];
Fig1 = figure();set(gcf,'Position',[450 50 1100 900]);
tiledlayout(4,4,"TileSpacing","compact","Padding","compact")

for i = 1:4
    ax = nexttile;
    imagesc(amide_map(:,:,i));axis off;colorbar;colormap(ax,colormap_Inferno_bbg)
    title('Amide '+tag(i));
    clim(amide_ran(i,:))

end

for i = 1:4
    ax = nexttile;
    imagesc(guan_map(:,:,i));axis off;colorbar;colormap(ax,colormap_Jet_black_bg)
    title('Guan '+tag(i));
    clim(guan_ran(i,:))

end

for i = 1:4
    ax = nexttile;
    imagesc(noe_map(:,:,i));axis off;colorbar;colormap(ax,colormap_Parula_black_bg)
    title('rNOE '+tag(i));
    clim(noe_ran(i,:))
end

nexttile;axis off
nexttile;axis off

for i = 3:4
    ax = nexttile;
    imagesc(mt_map(:,:,i-2));axis off;colorbar;colormap(ax,colormap_Jet_black_bg)
    title('MT '+tag(i));
    clim(mt_ran(i-2,:))
end