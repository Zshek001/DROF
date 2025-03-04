% display GE result of MTR,LDA,MPLF,DROF (single slice)
%   20241028
clear; clc
% close all; 
addpath ..\toolbox

load('20241029-colormap-black-bg.mat','colormap_Jet_black_bg','colormap_Parula_black_bg','colormap_Turbo_black_bg','colormap_Inferno_bbg')
%% load data
load('digital_brain.mat','roi_csf','roi_gm','roi_wm');

filepath = ".\20250204-noiseresult";
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
    load(filepath+"\"+filenameList(i));
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

%% regional analysis
index_csf = find(roi_csf==1); ncsf = length(index_csf);
index_gm = find(roi_gm==1); ngm = length(index_gm);
index_wm = find(roi_wm==1); nwm = length(index_wm);
amide_vec = reshape(amide_map,[],4);
guan_vec = reshape(guan_map,[],4);
rNOE_vec = reshape(noe_map,[],4);
MT_vec = reshape(mt_map,[],2);

amide_cell = cell(4,3); % three regions for four methods
guan_cell = cell(4,3);
rNOE_cell = cell(4,3);
MT_cell = cell(2,3);
for i = 1:4 % [MTR, LDA, MPLF, DROF]
    amide_cell{i,1} = amide_vec(index_csf,i);
    amide_cell{i,2} = amide_vec(index_gm,i);
    amide_cell{i,3} = amide_vec(index_wm,i);

    guan_cell{i,1} = guan_vec(index_csf,i);
    guan_cell{i,2} = guan_vec(index_gm,i);
    guan_cell{i,3} = guan_vec(index_wm,i);

    rNOE_cell{i,1} = rNOE_vec(index_csf,i);
    rNOE_cell{i,2} = rNOE_vec(index_gm,i);
    rNOE_cell{i,3} = rNOE_vec(index_wm,i);

    if i >= 3
        MT_cell{i-2,1} = MT_vec(index_csf,i-2);
        MT_cell{i-2,2} = MT_vec(index_gm,i-2);
        MT_cell{i-2,3} = MT_vec(index_wm,i-2);
    end
end

%% display
tag = ["\Delta Z [%]","\Delta Z [%]","\Delta Z [%]","R_{ex} [ms^{-1}]"];
amide_ran = [[-6,2]; [0,25]; [0,7]; [0,150]];
guan_ran = [[-2,0.5]; [0,16]; [0,10]; [0,75]];
noe_ran = [[-1,6]; [0,23]; [0,13.5]; [0,230]];
mt_ran = [[0,21]; [0,310]];
Fig1 = figure();set(gcf,'Position',[150 50 1100 900]);
tiledlayout(4,4,"TileSpacing","compact","Padding","compact")


for idx_method = 1:4
    ax = nexttile;
    data = [amide_cell{idx_method,1}; amide_cell{idx_method,2}; amide_cell{idx_method,3}];
    group = [repmat({'CSF'}, ncsf, 1); repmat({'GM'}, ngm, 1); repmat({'WM'}, nwm, 1)];

    boxplot(data, group);
    
    % title('Amide');
    ylabel(tag(idx_method),'FontWeight','bold');
end

for idx_method = 1:4
    ax = nexttile;
    data = [guan_cell{idx_method,1}; guan_cell{idx_method,2}; guan_cell{idx_method,3}];
    group = [repmat({'CSF'}, ncsf, 1); repmat({'GM'}, ngm, 1); repmat({'WM'}, nwm, 1)];

    boxplot(data, group);
    
    % title('Guan');
    ylabel(tag(idx_method),'FontWeight','bold');
end

for idx_method = 1:4
    ax = nexttile;
    data = [rNOE_cell{idx_method,1}; rNOE_cell{idx_method,2}; rNOE_cell{idx_method,3}];
    group = [repmat({'CSF'}, ncsf, 1); repmat({'GM'}, ngm, 1); repmat({'WM'}, nwm, 1)];

    boxplot(data, group);
    
    % title('rNOE');
    ylabel(tag(idx_method),'FontWeight','bold');
end

nexttile;axis off
nexttile;axis off

for idx_method = 3:4
    ax = nexttile;
    data = [MT_cell{idx_method-2,1}; MT_cell{idx_method-2,2}; MT_cell{idx_method-2,3}];
    group = [repmat({'CSF'}, ncsf, 1); repmat({'GM'}, ngm, 1); repmat({'WM'}, nwm, 1)];

    boxplot(data, group);
    
    % title('MT');
    ylabel(tag(idx_method),'FontWeight','bold');
end



%% save
resol = '600';
OutputDir = ".\Out";
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
fprintf("save to folder: "+join( split(OutputDir,'\'),'\\')+"\n");
exportgraphics(Fig1, OutputDir+"fig8.png", 'BackgroundColor', 'white', 'Resolution', resol);
