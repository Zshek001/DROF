% compare fitting results of best and worst fitting initial value
%   MPLF, MPLF2, PLOF, SROF, DROF
%   20250121 Huabin Zhang, huabinz@connect.hku.hk
%       compare all Delta Z value
%       only scatter plot
clear
% close all

%% calculate contribution and sensitivity
list = ["output_CO_MPLF_best.mat","output_CO_MPLF_2step_best.mat","output_DZ_PLOF_best.mat","output_DZ_SROF_best.mat","output_DZ_DROF_best.mat"....
        "output_CO_MPLF_worst.mat","output_CO_MPLF_2step_worst.mat","output_DZ_PLOF_worst.mat","output_DZ_SROF_worst.mat","output_DZ_DROF_worst.mat"];
folder = ".\BesAWor-ran5000-MS50\";
tag = {'MPLF_b','MPLF2_b','PLOF_b','SROF_b','DROF_b'...
         'MPLF_w','MPLF2_w','PLOF_w','SROF_w','DROF_w'   };

%% scatter plot
Fig1 = figure();set(gcf,'Position',[150 150 1000 450]);
tiledlayout(2,5,"TileSpacing","compact","Padding","compact")
for i = 1:length(list)
    nexttile
    load(folder+list(i),'model_4var_amide');
    mdl = model_4var_amide;

    legend_on = (i==1);
    titletag = tag(i);
    xtag = 'f: Amide [mM]';
    if any(i == [3,4,5,8,9,10], "all")
        ytag = 'R_{ex} [ms^{-1}]';
    else
        ytag = '\Delta Z [%]';
    end

    linearRegScatterPlot(mdl.Variables{:,2},mdl.Variables{:,end},xtag,ytag,titletag,legend_on);
end

Fig2 = figure();set(gcf,'Position',[150 150 1000 450]);
tiledlayout(2,5,"TileSpacing","tight","Padding","compact")
for i = 1:length(list)
    nexttile
    load(folder+list(i),'model_4var_guan');
    mdl = model_4var_guan;

    legend_on = (i==1);
    titletag = tag(i);
    xtag = 'f: CEST@2ppm [mM]';
    if any(i == [3,4,5,8,9,10], "all")
        ytag = 'R_{ex} [ms^{-1}]';
    else
        ytag = '\Delta Z [%]';
    end

    linearRegScatterPlot(mdl.Variables{:,1},mdl.Variables{:,end},xtag,ytag,titletag,legend_on);
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
exportgraphics(Fig1, OutputDir+"fig2_1.png", 'BackgroundColor', 'white', 'Resolution', resol);
exportgraphics(Fig2, OutputDir+"fig2_2.png", 'BackgroundColor', 'white', 'Resolution', resol);


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
    ylim([min(Y),1.2*max(Y)])

    p = plot(X1_sorted, Y_fitted_sorted, 'r-', 'LineWidth', 2);
    xlabel(xtag,'FontWeight','bold'); ylabel(ytag,'FontWeight','bold'); 
    % title(titletag);
    text(0.05, 0.93, ['R^2 = ' num2str(R2, '%.3f')], 'Units', 'normalized', 'Color', 'red', 'FontSize', 10,'FontWeight','bold');
    if legend_on
        lgd = legend('data', 'Fitted','Location','northwest');
        lgd.Position(1) = lgd.Position(1) - 0.006;
        lgd.Position(2) = lgd.Position(2) - 0.03; 
    end
    hold off;
end



