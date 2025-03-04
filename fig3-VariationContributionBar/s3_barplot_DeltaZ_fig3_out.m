% plot barchart based on Delta Z output
%   20250121 Huabin Zhang, huabinz@connect.hku.hk
%   no scatter plot
clear
% close all

%% calculate contribution and sensitivity
% list = ["output_ref.mat","output_DZ_MTR.mat","output_DZ_LDA.mat","output_DZ_MPLF.mat","output_DZ_PLOF.mat","output_DZ_MPLF_ss_2step.mat"];
list = ["output_DZ_MTR.mat","output_DZ_LDA.mat","output_DZ_LDA_2pool.mat","output_CO_MPLF.mat","output_CO_MPLF_2step.mat",...
    "output_R1rho_PLOF.mat","output_R1rho_SROF.mat","output_R1rho_DROF.mat"];
folder = ".\ran5000-fig3\";
tag = {'MTR_{asym}','LDA','LDA2','MPLF','DMPLF','PLOF','SROF','DROF'};

contri_amide_all = [];
contri_guan_all = [];
sens_amide_all = zeros(length(list),2); % slope and std
sens_guan_all = zeros(length(list),2);
sens_amide_R1rho = zeros(length(list),2); % for R1rho
sens_guan_R1rho = zeros(length(list),2);
for idx = 1:length(list)
    load(folder+list(idx),'model_4var_amide','model_4var_guan');
    
    % contribution based on ANOVA
    anovaTable = anova(model_4var_amide);
    contri_amide_all = [contri_amide_all, 100 * anovaTable.SumSq(1:4) / sum(anovaTable.SumSq)]; % one column for one algorithm

    anovaTable = anova(model_4var_guan);
    contri_guan_all = [contri_guan_all, 100 * anovaTable.SumSq(1:4) / sum(anovaTable.SumSq)]; % one column for one algorithm

    % sensitivity based on linear regression slope, unit [%/mM]
    if idx <= 5
        % output is Z-spec unit
        DeltaZ_guan = model_4var_guan.Variables{:,end};
        DeltaZ_amide = model_4var_amide.Variables{:,end};
    else
        % output is R1rho, should convert to Z-spec unit
        load(folder+list(idx),'Zamide_DZ_vec','Zguan_DZ_vec');
        DeltaZ_guan = Zguan_DZ_vec;
        DeltaZ_amide = Zamide_DZ_vec;

        % R1rho
        R1rho_guan = model_4var_guan.Variables{:,end};
        R1rho_amide = model_4var_amide.Variables{:,end};

        linFit = fitlm(model_4var_guan.Variables{:,1},R1rho_guan);
        sens_guan_R1rho(idx,1) = linFit.Coefficients.Estimate(2); % [ms^{-1}/mM]
        sens_guan_R1rho(idx,2) = linFit.Coefficients.SE(2);
    
        linFit = fitlm(model_4var_amide.Variables{:,2},R1rho_amide);
        sens_amide_R1rho(idx,1) = linFit.Coefficients.Estimate(2);
        sens_amide_R1rho(idx,2) = linFit.Coefficients.SE(2);

    end
    linFit = fitlm(model_4var_guan.Variables{:,1},DeltaZ_guan);
    sens_guan_all(idx,1) = linFit.Coefficients.Estimate(2); % [%/mM]
    sens_guan_all(idx,2) = linFit.Coefficients.SE(2);

    linFit = fitlm(model_4var_amide.Variables{:,2},DeltaZ_amide);
    sens_amide_all(idx,1) = linFit.Coefficients.Estimate(2);
    sens_amide_all(idx,2) = linFit.Coefficients.SE(2);
    
end

%% draw contribution and sensitivity (Delta Z)
barthick = 0.6;
legendfontsize = 8;
colorset=[[085,059,148];[152,114,202];[246,231,237];[216,163,152]];
ylim_amide = [0,140];
ylim_guan = [0,140];
ylim_amide_sen = [0.01,0.03]; % [%/mM]
ylim_guan_sen = [0.01,0.03]; % [%/mM]
indexshow = 1:8;

Fig1 = figure();set(gcf,'Position',[150 150 1000 450]);
tiledlayout(1,2,"TileSpacing","compact","Padding","loose")
nexttile

% 1-1 amide contribution
contri_amide_show = contri_amide_all([2,1,3,4],:); % resort display order
colorset_amide = colorset([2,1,3,4],:);

yyaxis left;
a = bar(indexshow,contri_amide_show(:,indexshow)',barthick, 'stacked'); 
ylim(ylim_amide)
% title('Variable Contribution and Sensitivity on Amide');
set(gca, 'xticklabel', tag(indexshow),'FontWeight','bold');
ylabel('Variance Contribution [%]','FontSize',15);

for i = 1:length(a)
    set(a(i),'FaceColor',colorset_amide(i,:)/255);
end

% 1-2 amide sensitivity
hold on;
yyaxis right;
plot(indexshow, sens_amide_all(indexshow,1), '-o', 'LineWidth', 2);
ylim(ylim_amide_sen);

lgd = legend({'Amide', 'CEST@2ppm', 'MT', 'rNOE', 'Sensitivity'}, 'Location', 'northeast','Orientation','vertical','FontSize',legendfontsize); 

% 2-1 guan contribution
nexttile
yyaxis left;
a = bar(indexshow,contri_guan_all(:,indexshow)',barthick, 'stacked'); 
ylim(ylim_guan)
% title('Variable Contribution and Sensitivity on Guan');
set(gca, 'xticklabel', tag(indexshow),'FontWeight','bold');

for i = 1:length(a)
    set(a(i),'FaceColor',colorset(i,:)/255);
end

% 2-2 guan sensitivity
hold on;
yyaxis right;
plot(indexshow, sens_guan_all(indexshow,1), '-o', 'LineWidth', 2);
ylabel('Sensitivity [%/mM]','FontSize',15);
ylim(ylim_guan_sen);

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
exportgraphics(Fig1, OutputDir+"fig3_1.png", 'BackgroundColor', 'white', 'Resolution', resol);
