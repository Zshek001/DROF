
% close all

% ----- Cross-platform path handling -----
here = fileparts(mfilename('fullpath'));  % folder of this script
if isempty(here)                          % e.g., if run as a Live Script
    here = pwd;
end

dataFolder = fullfile(here, '.', 'nB1-nTs_7');

%% load data
loadPath = fullfile(dataFolder, 'B1_Ts_List.mat');
load(loadPath,'B1List','TsList');
% B1List(end-2:end) = [];
nB1 = length(B1List);
nTs = length(TsList);

% sort file
filelListPath = fullfile(dataFolder, 'simData_B1_*.mat');
fileList = dir(filelListPath);
numValues = zeros(size(fileList));
for i = 1:length(fileList)
    tokens = regexp(fileList(i).name, '_(\d+\.?\d*)\.mat$', 'tokens');
    numValues(i) = str2double(tokens{1}{1});
end
[~, sortIdx] = sort(numValues);
fileList = fileList(sortIdx);

% initialize matrix
load(fullfile(fileList(1).folder, fileList(1).name),'sensR2_amide_all_givenB1');
[nx, ny, ~] = size(sensR2_amide_all_givenB1);

sensArr_all = zeros(3, 6, nTs, nB1); % 3 pools, 6 algorithms
R2Arr_all = zeros(3, 6, nTs, nB1); 
diffArr_all = zeros(3, 6, nTs, nB1); % 3pools, 6 algo

% load all data
for i = 1:nB1
    data = load(fullfile(fileList(i).folder, fileList(i).name),...
        'sensR2_amide_all_givenB1','sensR2_guan_all_givenB1','sensR2_noe_all_givenB1',...
        'diffss_meanList_all_givenB1');
        % sensR2_[pool]_all_givenB1: [2, 6, nTs], 2 for sensitivity and R2, 6 for algorithms 
    sensArr_all(1,:,:,i) = data.sensR2_amide_all_givenB1(1,:,:);
    sensArr_all(2,:,:,i) = data.sensR2_guan_all_givenB1(1,:,:);
    sensArr_all(3,:,:,i) = data.sensR2_noe_all_givenB1(1,:,:);

    R2Arr_all(1,:,:,i) = data.sensR2_amide_all_givenB1(2,:,:);
    R2Arr_all(2,:,:,i) = data.sensR2_guan_all_givenB1(2,:,:);
    R2Arr_all(3,:,:,i) = data.sensR2_noe_all_givenB1(2,:,:);

    diffArr_all(:,:,:,i) = data.diffss_meanList_all_givenB1(1:3,:,:);
end

% rearrange
[TsList,idxTsList] = unique(TsList);
R2Arr_all = R2Arr_all(:,:,idxTsList,:);
sensArr_all = sensArr_all(:,:,idxTsList,:);
TsList(end) = [];
R2Arr_all(:,:,end,:) = [];
sensArr_all(:,:,end,:) = [];

%% (1) plot sensitivity
Pre = {'amide', 'CEST@2ppm', 'rNOE'}; 
Post = {'MTRasym', 'LDA', 'DMPLF', 'QMTRasym', 'QLDA', 'QDMPLF'}; 
x_ticks = B1List;
y_ticks = TsList;
x_ticks_show = [0.4,0.8,1.2];
y_ticks_show = [1.0,1.5,2.0,2.5];

data = permute(sensArr_all,[4,3,1,2]); % [nB1, nTs, npool, nmethods]

Fig1 = figure();set(gcf,'Position',[50 50 1200 850]);
tiledlayout(3,6,"TileSpacing","loose","Padding","loose")
colormap parula; 
clim_pool = [[0,0.025];...
             [0,0.025];...
             [0,0.006]];

for idxpool = 1:3
    for idxmeth = 1:6
        nexttile
        imagesc(x_ticks, y_ticks, squeeze(data(:,:,idxpool,idxmeth))'); 
        clim(clim_pool(idxpool,:))
        set(gca, 'YDir', 'normal','FontWeight','bold','FontSize',12,'TickDir','out')
        set(gca, 'XTick', x_ticks_show,'FontWeight','normal','FontSize',10);
        set(gca, 'YTick', y_ticks_show,'FontWeight','normal','FontSize',10);

        % for output, discard text
        xlabel('B1 [uT]'); ylabel('Ts [s]');
        title(sprintf('%s - %s', Pre{idxpool}, Post{idxmeth}));
    end
end

% create colorbar for each row
ypos = [0.71,0.415,0.12];
for idxpool = 1:3
    ax = nexttile(idxpool*6);
    
    cbh = colorbar(ax);
    cbh.Layout.Tile = 'east';
    
    pos = cbh.Position;
    pos(2) = ypos(idxpool); % y position
    pos(3) = 0.02; % width
    pos(4) = 0.2; % height
    set(cbh,'Location','manual','Position', pos); % width can be changed only in 'manual' mode
    cbh.FontSize = 12; 
    cbh.FontWeight = 'bold'; 
end

%% (2) plot R2
data = permute(R2Arr_all,[4,3,1,2]); % [nB1, nTs, npool, nmethods]

Fig2 = figure();set(gcf,'Position',[50 50 1200 850]);
tiledlayout(3,6,"TileSpacing","loose","Padding","loose")
colormap jet; 

for idxpool = 1:3
    for idxmeth = 1:6
        nexttile
        imagesc(x_ticks, y_ticks, squeeze(data(:,:,idxpool,idxmeth))'); 
        clim([0,1])
        set(gca, 'YDir', 'normal','FontWeight','bold','FontSize',12,'TickDir','out')
        set(gca, 'XTick', x_ticks_show,'FontWeight','normal','FontSize',10);
        set(gca, 'YTick', y_ticks_show,'FontWeight','normal','FontSize',10);

        % for output, discard text
        xlabel('B1 [uT]'); ylabel('Ts [s]');
        title(sprintf('%s - %s', Pre{idxpool}, Post{idxmeth}));
    end
end

cbh = colorbar; 
cbh.Layout.Tile = 'east';

temp = cbh.Position;
temp(3) = 0.03;
set(cbh,'Location','manual','Position',temp) % width can be changed only in 'manual' mode

cbh.FontSize = 12; 
cbh.FontWeight = 'bold'; 

%% (3) plot R2
data = permute(abs(diffArr_all),[4,3,1,2]); % [nB1, nTs, npool, nmethods]

Fig2 = figure();set(gcf,'Position',[50 50 1200 850]);
tiledlayout(3,6,"TileSpacing","loose","Padding","loose")
colormap jet; 

for idxpool = 1:3
    for idxmeth = 1:6
        nexttile
        imagesc(x_ticks, y_ticks, squeeze(data(:,:,idxpool,idxmeth))'); 
        % clim([0,1])
        set(gca, 'YDir', 'normal','FontWeight','bold','FontSize',12,'TickDir','out')
        set(gca, 'XTick', x_ticks_show,'FontWeight','normal','FontSize',10);
        set(gca, 'YTick', y_ticks_show,'FontWeight','normal','FontSize',10);

        % for output, discard text
        xlabel('B1 [uT]'); ylabel('Ts [s]');
        title(sprintf('%s - %s', Pre{idxpool}, Post{idxmeth}));
    end
end

cbh = colorbar; 
cbh.Layout.Tile = 'east';

temp = cbh.Position;
temp(3) = 0.03;
set(cbh,'Location','manual','Position',temp) % width can be changed only in 'manual' mode

cbh.FontSize = 12; 
cbh.FontWeight = 'bold'; 
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
% exportgraphics(Fig1, OutputDir+"fig_mB1_mTs_sen.png", 'BackgroundColor', 'white', 'Resolution', resol);
% exportgraphics(Fig2, OutputDir+"fig_mB1_mTs_R2.png", 'BackgroundColor', 'white', 'Resolution', resol);
