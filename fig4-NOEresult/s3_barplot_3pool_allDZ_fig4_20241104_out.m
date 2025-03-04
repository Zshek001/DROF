% plot barchart based on Delta Z output
%   20241024 Huabin Zhang, huabinz@connect.hku.hk
%       compare amide, guan, NOE
%       compare two offs
%       compare MTR, LDA, MPLF, MPLF2-ss
%   20241101: compare all delta Z, remove sensitivity
clear
% close all

%% calculate contribution and sensitivity
list = ["output_DZ_MTR.mat","output_DZ_LDA.mat","output_CO_MPLF.mat","output_R1rho_DROF.mat"];
folder = ".\ran5000-fig4\";
tag = {'MTR_{asym}','LDA','MPLF','DROF'};

contri_amide_all = [];
contri_guan_all = [];
contri_noe_all = [];
sens_amide_all = zeros(length(list),2); % slope and std
sens_guan_all = zeros(length(list),2);
sens_noe_all = zeros(length(list),2);
for idx = 1:length(list)
    load(folder+list(idx),'model_4var_amide','model_4var_guan','model_4var_noe');
    
    % contribution based on ANOVA
    anovaTable = anova(model_4var_amide);
    contri_amide_all = [contri_amide_all, 100 * anovaTable.SumSq(1:4) / sum(anovaTable.SumSq)]; % one column for one algorithm

    anovaTable = anova(model_4var_guan);
    contri_guan_all = [contri_guan_all, 100 * anovaTable.SumSq(1:4) / sum(anovaTable.SumSq)];

    anovaTable = anova(model_4var_noe);
    contri_noe_all = [contri_noe_all, 100 * anovaTable.SumSq(1:4) / sum(anovaTable.SumSq)];

    % sensitivity based on linear regression slope, unit [%/mM]
    if idx <= 3
        % output is Z-spec unit
        DeltaZ_guan = model_4var_guan.Variables{:,end};
        DeltaZ_amide = model_4var_amide.Variables{:,end};
        DeltaZ_noe = model_4var_noe.Variables{:,end};
    else
        % output is R1rho, should convert to Z-spec unit
        load(folder+list(idx),'Zamide_DZ_vec','Zguan_DZ_vec','Znoe_DZ_vec','fit_para');
        DeltaZ_guan = Zguan_DZ_vec;
        DeltaZ_amide = Zamide_DZ_vec;
        DeltaZ_noe = Znoe_DZ_vec;

    end
    linFit = fitlm(model_4var_guan.Variables{:,1},DeltaZ_guan);
    sens_guan_all(idx,1) = linFit.Coefficients.Estimate(2); % [%/mM]
    sens_guan_all(idx,2) = linFit.Coefficients.SE(2);

    linFit = fitlm(model_4var_amide.Variables{:,2},DeltaZ_amide);
    sens_amide_all(idx,1) = linFit.Coefficients.Estimate(2);
    sens_amide_all(idx,2) = linFit.Coefficients.SE(2);

    linFit = fitlm(model_4var_noe.Variables{:,4},DeltaZ_noe);
    sens_noe_all(idx,1) = linFit.Coefficients.Estimate(2);
    sens_noe_all(idx,2) = linFit.Coefficients.SE(2);
    
end

%% draw contribution and sensitivity
barthick = 0.6;
legendfontsize = 8;
colorset=[[085,059,148];[152,114,202];[246,231,237];[216,163,152]];
ylim_amide = [0,140];
ylim_guan = [0,140];
ylim_noe = [0,140];
ylim_amide_sen = [0.01,0.03]; % [%/mM]
ylim_guan_sen = [0.01,0.03]; % [%/mM]
ylim_noe_sen = [0.01,0.03]; % [%/mM]

Fig2 = figure();set(gcf,'Position',[150 150 1000 450]);
tiledlayout(1,3,"TileSpacing","loose","Padding","loose")
nexttile

% 1-1 amide contribution
contri_amide_show = contri_amide_all([2,1,3,4],:); % resort display order
colorset_amide = colorset([2,1,3,4],:);

yyaxis left;
a = bar(1:length(list),contri_amide_show',barthick, 'stacked'); 
ylim(ylim_amide)
% title('Amide @ 3.5 ppm');
set(gca, 'xticklabel', tag,'FontSize',10,'FontWeight','bold');
ylabel('Variance Contribution [%]','FontSize',15);

for i = 1:length(a)
%     set(a(i),'FaceColor',colorset_amide(i,:)/255,'edgecolor','none');
    set(a(i),'FaceColor',colorset_amide(i,:)/255);
end

% 1-2 amide sensitivity
hold on;
yyaxis right;
plot(1:length(list), sens_amide_all(:,1), '-o', 'LineWidth', 2);
ylim(ylim_amide_sen);

lgd = legend({'Amide', 'CEST@2ppm', 'MT', 'rNOE', 'Sensitivity'}, 'Location', 'northeast','Orientation','vertical','FontSize',legendfontsize); 

% 2-1 guan contribution
nexttile
yyaxis left;
a = bar(1:length(list),contri_guan_all',barthick, 'stacked'); 
ylim(ylim_guan)
% title('Guan @ 2.0 ppm');
set(gca, 'xticklabel', tag,'FontSize',10,'FontWeight','bold');

for i = 1:length(a)
    set(a(i),'FaceColor',colorset(i,:)/255);
end

% 2-2 guan sensitivity
hold on;
yyaxis right;
plot(1:length(list), sens_guan_all(:,1), '-o', 'LineWidth', 2);
ylim(ylim_guan_sen);

% 3-1 noe contribution
nexttile
contri_noe_show = contri_noe_all([4,2,3,1],:); % resort display order
colorset_noe = colorset([4,2,3,1],:);

yyaxis left;
a = bar(1:length(list),contri_noe_show',barthick, 'stacked'); 
ylim(ylim_noe)
% title('NOE @ -3.5ppm');
set(gca, 'xticklabel', tag,'FontSize',10,'FontWeight','bold');

for i = 1:length(a)
    set(a(i),'FaceColor',colorset_noe(i,:)/255);
end

% 3-2 noe sensitivity
hold on;
yyaxis right;
plot(1:length(list), sens_amide_all(:,1), '-o', 'LineWidth', 2);
ylabel('Sensitivity [%/mM]','FontSize',15);
ylim(ylim_noe_sen);


%% scatter plot
Fig1 = figure();set(gcf,'Position',[150 150 1000 450]);
tiledlayout('flow',"TileSpacing","compact","Padding","compact")
for i = 1:length(list)
    nexttile
    load(folder+list(i),'model_4var_amide');
    mdl = model_4var_amide;

    legend_on = 0;
    titletag = tag(i);
    xtag = 'f: Amide [mM]';
    if i == 4
        ytag = 'R_{ex} [ms^{-1}]';
    else
        ytag = '\Delta Z [%]';
    end

    linearRegScatterPlot(mdl.Variables{:,2},mdl.Variables{:,end},xtag,ytag,titletag,legend_on);
end

for i = 1:length(list)
    nexttile
    load(folder+list(i),'model_4var_guan');
    mdl = model_4var_guan;

    legend_on = 0;
    titletag = tag(i);
    xtag = 'f: CEST@2ppm [mM]';
    if i == 4
        ytag = 'R_{ex} [ms^{-1}]';
    else
        ytag = '\Delta Z [%]';
    end

    linearRegScatterPlot(mdl.Variables{:,1},mdl.Variables{:,end},xtag,ytag,titletag,legend_on);
end

for i = 1:length(list)
    nexttile
    load(folder+list(i),'model_4var_noe');
    mdl = model_4var_noe;

    legend_on = (i==4);
    titletag = tag(i);
    xtag = 'f: rNOE [mM]';
    if i == 4
        ytag = 'R_{ex} [ms^{-1}]';
    else
        ytag = '\Delta Z [%]';
    end

    linearRegScatterPlot(mdl.Variables{:,4},mdl.Variables{:,end},xtag,ytag,titletag,legend_on);
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
exportgraphics(Fig1, OutputDir+"fig4_1.png", 'BackgroundColor', 'white', 'Resolution', resol);
exportgraphics(Fig2, OutputDir+"fig4_2.png", 'BackgroundColor', 'white', 'Resolution', resol);

%%
function p = linearRegScatterPlot(X,Y,xtag,ytag,titletag,legend_on)
    % linear regression
    linFit = fitlm(X, Y);
    Y_fitted = linFit.Fitted;
    residuals = Y - Y_fitted;
    
    % 95% confidence boundary
    n = length(X);
    SE = sqrt(sum(residuals.^2) / (n - 2)); 
    t_val = tinv(0.975, n - 2); 
    conf_interval = t_val * SE * sqrt(1/n + (X - mean(X)).^2 / sum((X - mean(X)).^2));
    R2 = linFit.Rsquared.Ordinary;
    
    % sort data for plotting
    [X1_sorted, idx] = sort(X);
    Y_fitted_sorted = Y_fitted(idx);
    conf_interval_sorted = conf_interval(idx);

    % scatter plot without CB
    scatter(X, Y, 2, 'b', 'filled'); hold on;
    p = plot(X1_sorted, Y_fitted_sorted, 'r-', 'LineWidth', 2);
    xlabel(xtag,'FontWeight','bold'); ylabel(ytag,'FontWeight','bold'); 
    % title(titletag);
    text(0.05, 0.93, ['R^2 = ' num2str(R2, '%.3f')], 'Units', 'normalized', 'Color', 'red', 'FontSize', 10,'FontWeight','bold');
    if legend_on
        lgd = legend('data', 'Fitted','Location','southeast');
        lgd.Position(1) = lgd.Position(1) + 0.01;
        lgd.Position(2) = lgd.Position(2) - 0.02;
    end
    hold off;
end



