% display GE result of MTR,LDA,MPLF,DROF (single slice)
%   20241028
clear; clc
% close all; 
addpath ..\toolbox

load('20241029-colormap-black-bg.mat','colormap_Jet_black_bg','colormap_Parula_black_bg','colormap_Turbo_black_bg','colormap_Inferno_bbg')
%% load data
filepath = uigetdir(".\","select CEST fitting results folder");
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

%% display
tag = ["Z [%]","Z [%]","Z [%]","R1rho [ms^{-1}]"];
amide_ran = [[-6,2]; [0,25]; [0,7]; [0,150]];
guan_ran = [[-2,0.5]; [0,16]; [0,10]; [0,75]];
noe_ran = [[-1,6]; [0,23]; [0,13.5]; [0,230]];
mt_ran = [[0,21]; [0,310]];
Fig1 = figure();set(gcf,'Position',[150 50 1100 900]);
tiledlayout(4,4,"TileSpacing","compact","Padding","compact")


for i = 1:4
    ax = nexttile;
    imagesc(amide_map(:,:,i));axis off;colorbar;colormap(ax,colormap_Inferno_bbg)
    % title('Amide '+tag(i));
    clim(amide_ran(i,:))

end

for i = 1:4
    ax = nexttile;
    imagesc(guan_map(:,:,i));axis off;colorbar;colormap(ax,colormap_Jet_black_bg)
    % title('Guan '+tag(i));
    clim(guan_ran(i,:))

end


for i = 1:4
    ax = nexttile;
    imagesc(noe_map(:,:,i));axis off;colorbar;colormap(ax,colormap_Parula_black_bg)
    % title('rNOE '+tag(i));
    clim(noe_ran(i,:))
end

nexttile;axis off
nexttile;axis off

for i = 3:4
    ax = nexttile;
    imagesc(mt_map(:,:,i-2));axis off;colorbar;colormap(ax,colormap_Jet_black_bg)
    % title('MT '+tag(i));
    clim(mt_ran(i-2,:))
end

%% save
% resol = '600';
% OutputDir = ".\Out";
% DirTemp = OutputDir;
% cnt = 1;
% while(1)
%     if exist(DirTemp,'dir')
%         DirTemp = OutputDir + "_" + num2str(cnt);
%         cnt = cnt + 1;
%     else
%         break;
%     end
% end
% OutputDir = DirTemp + "\";
% mkdir(OutputDir);
% fprintf("save to folder: "+join( split(OutputDir,'\'),'\\')+"\n");
% exportgraphics(Fig1, OutputDir+"fig7.png", 'BackgroundColor', 'white', 'Resolution', resol);
