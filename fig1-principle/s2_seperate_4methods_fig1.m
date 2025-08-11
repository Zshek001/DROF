% perform PLOF on phantom data (random setting)
% use multi-linear regression to assess dependence on MTf, NOEf
% 2024-1101
clear
warning off

% ----- Cross-platform path handling -----
here = fileparts(mfilename('fullpath'));  % folder of this script
if isempty(here)                          % e.g., if run as a Live Script
    here = pwd;
end

% Add toolbox folder (.. / toolbox)
toolboxDir = fullfile(here, '..', 'toolbox');
if exist(toolboxDir, 'dir')
    addpath(toolboxDir);
else
    warning('Toolbox folder not found: %s', toolboxDir);
end

% If you need myfcn subfolder inside toolbox, uncomment:
% myfcnDir = fullfile(here, 'toolbox', 'myfcn');
% if exist(myfcnDir, 'dir'); addpath(myfcnDir); end

display_on = 1;
multistart_N = 20; % times for MultiStart

% Data file alongside this script
zspec_file = fullfile(here, 'zspec_user1demo.mat');

%% load data
load(char(zspec_file), "offs", "zspec"); % m0 not included
    % offs: [nf,1]
    % zspec: [nf, 1]
[nf,~] = size(zspec);

Fig1 = figure(); set(gcf,'Position',[150 50 1200 850]);
tiledlayout(2,2,"TileSpacing","loose","Padding","loose")
yran = [0,1];
xran = [-8,8];
xran_plof = [0,8];
legfont = 8;
maksiz = 6;
Linewd = 1.5;

%% (1) LDA
nexttile
[Z_ref, LD] = fitting_LDA_demo(offs, zspec, multistart_N);
indbg = find((offs<=0.8 & offs>=-0.8)|abs(offs)>=7.5);
indother = find((offs>0.8 | offs<-0.8) & abs(offs)<7.5);
hold on
plot(offs(indbg),zspec(indbg),'go','MarkerSize',maksiz,'LineWidth',Linewd);
plot(offs(indother),zspec(indother),'bo','MarkerSize',maksiz);
plot(offs,Z_ref,'b',offs,LD,'r-.','LineWidth',Linewd);
xlabel('offs [ppm]');ylabel('Z-value')
set(gca,'XDir','reverse', 'FontWeight', 'bold', 'FontSize', 14)
lgd = legend('Data (for fitting)','Data (for subtracting)','Z_{fit} (DS)','Diff.','location','southwest','fontsize',legfont);
lgd.Position(2) = lgd.Position(2) + 0.02;
ylim(yran); xlim(xran);
hold off

%% (2) MPLF
nexttile
[Z_fit,Z_DS, Z_amide, Z_MT, Z_noe, Z_guan] = fitting_MPLF_demo(offs, zspec, multistart_N);
hold on
plot(offs,zspec,'bo','MarkerSize',maksiz);
plot(offs,Z_fit,'b',offs,Z_DS,'r-.',offs,Z_amide,'g-.',offs,Z_MT,'m-.', ...
     offs,Z_noe,'c-.',offs,Z_guan,'k-.','LineWidth',Linewd)
hold off
xlabel('offs [ppm]');ylabel('Z-value')
set(gca,'XDir','reverse', 'FontWeight', 'bold', 'FontSize', 14)
legend('Data', 'Z_{fit} (5 pools)', 'Water', 'Amide', 'MT','rNOE','CEST@2ppm', 'Location', 'southeast', ...
    'fontsize',legfont);
ylim(yran); xlim(xran);

%% (3) PLOF
nexttile
B0 = 3; % T
FitParam.R1 = 1; % brain ~ 1s @ 3T
FitParam.satpwr = 0.8; % [uT]
FitParam.tsat = 2; % [s]
FitParam.Magfield = 42.5764 * B0; % [Hz]
FitParam.PeakOffset = 3.5; % used in PLOF, position to be polynomialize

indroi = find( offs <= 8 & offs >= 0.5 );
offs_plof = offs(indroi);
indbg = find( offs_plof <= 1.5 | offs_plof >= 5 );
indcest = find( offs_plof > 1.5 & offs_plof < 5 );
[Z_fit, Z_bak, Z_amide, Z_guan] = fitting_PLOF_demo(offs_plof, zspec, FitParam, multistart_N);
hold on
plot(offs(indbg),zspec(indbg),'go','MarkerSize',maksiz,'LineWidth',Linewd);
plot(offs(indcest),zspec(indcest),'bo','MarkerSize',maksiz);
plot(offs_plof,Z_bak,'b',offs_plof,Z_fit,'k',offs_plof,Z_amide,'m:', ...
     offs_plof,Z_guan,'k:','LineWidth',Linewd)
hold off
xlabel('offs [ppm]');ylabel('Z-value')
set(gca,'XDir','reverse', 'FontWeight', 'bold', 'FontSize', 14)
legend('Data (for 1st fitting)','Data (for 2nd fitting)', 'Z_{fit1} (DS+MT)','Z_{fit2} (all)', 'Amide','CEST@2ppm', 'Location', 'southwest', ...
    'fontsize',legfont);
ylim(yran); xlim(xran_plof);

%% (4) DROF
nexttile
indbg = find( offs <= 1 | offs >= 6 );
indcest = find( offs > 1 & offs < 6 );

%[Z_fit1, Z_fit2, Z_DS, Z_amide, Z_MT, Z_noe, Z_guan] = fitting_MPLF_2step_demo(offs, zspec, multistart_N);
[Z_fit1, Z_fit2, Z_bak, Z_guan, Z_amide] = fitting_DROF_demo(offs, zspec, FitParam, multistart_N);

hold on
plot(offs(indbg),zspec(indbg),'go','MarkerSize',maksiz,'LineWidth',Linewd);
plot(offs(indcest),zspec(indcest),'bo','MarkerSize',maksiz);
plot(offs,Z_fit1,'b',offs,Z_fit2,'k',offs,Z_bak,'r:',offs,Z_amide,'m:', ...
     offs,Z_guan,'k:','LineWidth',Linewd)
hold off
xlabel('offs [ppm]');ylabel('Z-value')
set(gca,'XDir','reverse', 'FontWeight', 'bold', 'FontSize', 14)
legend('Data (for 1st fitting)','Data (for 2nd fitting)', 'Z_{fit1} (DS+MT+rNOE)', 'Z_{fit2} (5 pools)', 'Water+MT+NOE (background)', 'background + Amide', 'MT','background + CEST@2ppm', 'Location', 'southeast', ...
    'fontsize',legfont);
ylim(yran); xlim(xran);

%% save
resol = 600;
outBase = fullfile(here, 'Out');  % base output folder next to this script
outDir = outBase;
cnt = 1;
while exist(outDir, 'dir')
    outDir = sprintf('%s_%d', outBase, cnt);
    cnt = cnt + 1;
end
mkdir(outDir);
fprintf('save to folder: %s\n', outDir);
exportgraphics(Fig1, fullfile(outDir, 'fig1_1.png'), 'BackgroundColor', 'white', 'Resolution', resol);

%% fitting function
function [Z_ref, LD] = fitting_LDA_demo(offs, zspec_vec, multistart_N)

    %% LDA setting
    [nf,Npixel] = size(zspec_vec);
    fit_rgc = [-0.8, 0.8]; % Center range of Z-spectrum used for Lorentzian fitting, in ppm (0.5~1.5 ppm)
    fit_rge = [-19, 19]; % Edge range of Z-spectrum used for Lorentzian fitting, in ppm 

    %              1. water                    2. MT
    %      Zi   A1     G1      dw1       A2      G2      dw2   
    iv = [ 1    0.5    2       0        0.1      60      -2.5];
    lb = [ 0.5  0.02   0.3    -0.2      0.001    30      -3.0];
    ub = [ 1    1      100     +0.2      0.4      100     -2.0];

    iv = iv(1:4);
    lb = lb(1:4);
    ub = ub(1:4);

    fit_para = zeros(Npixel, length(iv));
    Zamide_DZ_vec = zeros(1,Npixel);
    Zguan_DZ_vec = zeros(1,Npixel);

    %% LDA fitting
    backNum = 0;
    tic
    for i = 1:Npixel
        fprintf(repmat('\b',1,backNum));
        backNum = fprintf('Calculation process: %d/%d',i,Npixel);

        zspec = zspec_vec(:,i);

        fit_ind = get_fit_ind(zspec, offs, fit_rgc, fit_rge);
        offs_trunc = offs(fit_ind);
        zspec_trunc = zspec(fit_ind);

        problem = createOptimProblem('lsqcurvefit','x0',iv,'objective',@lorentzMultipool,...
            'lb',lb,'ub',ub,'xdata',offs_trunc,'ydata',zspec_trunc);
        
        ms = MultiStart('Display','off');
        [par,~] = run(ms,problem,multistart_N);

        % fitted curves
        Z_ref = lorentzMultipool(par,offs);
        LD = Z_ref - zspec;
        [~, ind_amide] = min(abs(offs - 3.5));
        [~, ind_guan] = min(abs(offs - 2));

        Zamide_DZ_vec(i) = 100*LD(ind_amide);
        Zguan_DZ_vec(i) = 100*LD(ind_guan);
    end
    fprintf(repmat('\b',1,backNum));
    fprintf('elapsed time: %.4f s\n', toc);

end

function [Z_fit, Z_DS, Z_amide, Z_MT, Z_noe, Z_guan] = fitting_MPLF_demo(offs, zspec_vec, multistart_N)
% OUTPUT
%   Zamide_DZ_vec: [1,Npixel]
%   Zguan_DZ_vec: [1,Npixel]
%   fit_para: [16, Npixel]
%       Zi (1), Water (2-4), Amide (5-7), NOE (8-10), MT (11-13), Guan (14-16)

    %% MPLF setting
    [nf,Npixel] = size(zspec_vec);

    %             1. Water              2. Amide               3. NOE                 4. MT                  5. Guan
    %      Zi     A1    G1    dw1       A2     G2    dw2       A3     G3    dw3       A4     G4    dw4       A5     G5    dw5
    iv = [ 0.9    0.9   2.3   0         0.01   2.0   3.5       0.1    4.0   -3.5      0.1    60    -2.5      0.01   2.0   2.0];
    lb = [ 0.5    0.4   1.0   -0.2      0.0    0.5   3.4       0.0    1.0   -3.7      0.001  30    -3.0      0.0    0.5   1.8];
    ub = [ 1.0    1.0   6.0   +0.2      0.2    10.0  3.6       0.2    15.0  -3.3      0.4    100   -2.0      0.2    10.0  2.2];
    fit_para = zeros(length(iv), Npixel);
    Zamide_DZ_vec = zeros(1,Npixel); % based on delta Z
    Zguan_DZ_vec = zeros(1,Npixel);

    %% MPLF fitting
    backNum = 0;
    tic
    for i = 1:Npixel
        fprintf(repmat('\b',1,backNum));
        backNum = fprintf('Calculation process: %d/%d',i,Npixel);

        zspec = zspec_vec(:,i);

        problem = createOptimProblem('lsqcurvefit','x0',iv,'objective',@lorentzMultipool,...
            'lb',lb,'ub',ub,'xdata',offs,'ydata',zspec);
        ms = MultiStart('Display','off');
        [par,~] = run(ms,problem,multistart_N);

        % fitted curves
        Z_fit = lorentzMultipool(par,offs);
        Z_DS = lorentzMultipool([par(1:4)],offs); % Water
        Z_amide = lorentzMultipool([par(1),par(5:7)], offs); % amide
        Z_noe = lorentzMultipool([par(1),par(8:10)], offs); 
        Z_MT = lorentzMultipool([par(1),par(11:13)], offs); 
        Z_guan = lorentzMultipool([par(1),par(14:16)], offs); 

        Zamide_DZ_vec(i) = 100*par(5);
        Zguan_DZ_vec(i) = 100*par(14);
    end
    fprintf(repmat('\b',1,backNum));
    fprintf('elapsed time: %.4f s\n', toc);

end

function [Z_fit1, Z_fit2, Z_DS, Z_amide, Z_MT, Z_noe, Z_guan] = fitting_MPLF_2step_demo(offs, zspec_vec, multistart_N)
% INPUT
%   offs: [nf,1]
%   zspec_vec: [nf,Npixel]
% OUTPUT
%   Zamide_DZ_vec: [1,Npixel]
%   Zguan_DZ_vec: [1,Npixel]
%   fit_para: [16, Npixel]
%       Zi (1), Water (2-4), Amide (5-7), NOE (8-10), MT (11-13), Guan (14-16)

    %% MPLF setting
    [nf,Npixel] = size(zspec_vec);

    %             1. Water              2. Amide               3. NOE                 4. MT                  5. Guan
    %      Zi     A1    G1    dw1       A2     G2    dw2       A3     G3    dw3       A4     G4    dw4       A5     G5    dw5
    iv = [ 0.9    0.9   2.3   0         0.01   2.0   3.5       0.1    4.0   -3.5      0.1    60    -2.5      0.01   2.0   2.0];
    lb = [ 0.5    0.4   1.0   -1.0      0.0    0.5   3.1       0.0    1.0   -3.9      0.001  30    -3.0      0.0    0.5   1.5];
    ub = [ 1.0    1.0   6.0   +1.0      0.3    20.0  3.9       0.3    15.0  -3.1      0.4    100   -2.0      0.3    20.0  2.5];
    fit_para = zeros(length(iv), Npixel);
    Zamide_DZ_vec = zeros(1,Npixel); % based on delta Z
    Zguan_DZ_vec = zeros(1,Npixel);

    offs_bak_ind = find( offs <= 1 | offs >= 6 );
    offs_bak = offs(offs_bak_ind);
    fit_bak_ind = [1:4,8:13]; % Water + MT + NOE

    %% MPLF fitting
    backNum = 0;
    tic
    for i = 1:Npixel
        fprintf(repmat('\b',1,backNum));
        backNum = fprintf('Calculation process: %d/%d',i,Npixel);

        zspec = zspec_vec(:,i);

        %====================== first fitting (DS+MT) =======================
        problem = createOptimProblem('lsqcurvefit','x0',iv(fit_bak_ind),'objective',@lorentzMultipool,...
            'lb',lb(fit_bak_ind),'ub',ub(fit_bak_ind),'xdata',offs_bak,'ydata',zspec(offs_bak_ind));
        ms = MultiStart('Display','off');
        [par_bak,~] = run(ms,problem,multistart_N);

        %====================== second fitting (fix par_bak) =======================
        iv_temp = iv; iv_temp(fit_bak_ind) = par_bak;
        lb_temp = lb; lb_temp(fit_bak_ind) = par_bak;
        ub_temp = ub; ub_temp(fit_bak_ind) = par_bak;

        problem = createOptimProblem('lsqcurvefit','x0',iv_temp,'objective',@lorentzMultipool,...
            'lb',lb_temp,'ub',ub_temp,'xdata',offs,'ydata',zspec);
        ms = MultiStart('Display','off');
        [par,~] = run(ms,problem,multistart_N);
        fit_para(:,i) = par;

        % fitted curves
        Z_fit1 = lorentzMultipool(par_bak,offs);
        Z_fit2 = lorentzMultipool(par,offs);
        Z_DS = lorentzMultipool([par(1:4)],offs); % Water
        Z_amide = lorentzMultipool([par(1),par(5:7)], offs); % amide
        Z_noe = lorentzMultipool([par(1),par(8:10)], offs); 
        Z_MT = lorentzMultipool([par(1),par(11:13)], offs); 
        Z_guan = lorentzMultipool([par(1),par(14:16)], offs); 

        Zamide_DZ_vec(i) = 100*par(5);
        Zguan_DZ_vec(i) = 100*par(14);
    end
    fprintf(repmat('\b',1,backNum));
    fprintf('elapsed time: %.4f s\n', toc);

end

function [Z_fit, Z_bak, Z_amide, Z_guan] = fitting_PLOF_demo(offs, zspec_vec, FitParam, multistart_N)
% OUTPUT
%   fit_para: [8, Npixel]
%       Water (1-2), Amide (3-5), Guan (6-8)
    %% PLOF setting
    FitParam.PeakOffset = 3.5; % used in CurveFunction, position to be polynomialize
    FitParam.MultiStart = multistart_N;
    WholeRange = [0.5, 8];    % signal for background/CEST fitting
    PeakRange = [1.5, 5];     % extract background signal
    
    [nf,Npixel] = size(zspec_vec);
    Ramide_vec = zeros(1,Npixel);
    Rguan_vec = zeros(1,Npixel);
    
    Zamide_DZ_vec = zeros(1,Npixel);
    Zguan_DZ_vec = zeros(1,Npixel);
    Zamide_Co_vec = zeros(1,Npixel);
    Zguan_Co_vec = zeros(1,Npixel);
    fit_para = []; % [npool, Npixel]
    
    MTbg_vec = zeros(1,Npixel);
    
    %% PLOF fitting
    backNum = 0;
    tic
    for i = 1:Npixel
        fprintf(repmat('\b',1,backNum));
        backNum = fprintf('Calculation process: %d/%d',i,Npixel);
                  
        Z_spectrum = zspec_vec(:,i);
        Z_spectrum = squeeze(Z_spectrum);
    
        % only use zspec without m0
        [Z_fit, Z_bak, Z_amide, Z_guan] = PLOF_custom_demo(offs, Z_spectrum, WholeRange, PeakRange, FitParam);
    end
    fprintf(repmat('\b',1,backNum));
    fprintf('elapsed time: %.4f s\n', toc);
end

function [Z_fit1, Z_fit2, Z_bak, Z_guan, Z_amide] = fitting_DROF_demo(offs, zspec_vec, FitParam, multistart_N)
% use steady-state signal equation: "Z = (1 - cos2thet*R1/R1rho).* exp(-R1rho * tsat) + cos2thet * R1/R1rho"
% OUTPUT
%   Zamide_DZ_vec: [1,Npixel]
%   Zguan_DZ_vec: [1,Npixel]
%   fit_para: [16, Npixel]
%       Zi (1), Water (2-4), Amide (5-7), NOE (8-10), MT (11-13), Guan (14-16)

    %% MPLF setting
    [nf,Npixel] = size(zspec_vec);
    % % saturation settings
    % pulse1_pwr = 0.8; % in uT
    % pulse1_dur = 2; % pulse duration in s
    % gamma_hz = 42.5764;
    % B0 = 3;
    % FitParam.R1 = 1; % no T1 correction
    % FitParam.satpwr = pulse1_pwr; % saturation power (uT)
    % FitParam.tsat = pulse1_dur; % saturation length (second), 100s to make it steady-state
    % FitParam.Magfield = gamma_hz * B0;
    R1rho_ss = @(par,offs) lorentzMultipool_R1rho(par,offs,FitParam);

    %             1. Water              2. Amide               3. NOE                 4. MT                  5. Guan
    %      Zi     A1    G1    dw1       A2     G2    dw2       A3     G3    dw3       A4     G4    dw4       A5     G5    dw5
    iv = [ 0.9    3     1.0   0         0.1    0.5   3.5       0.1    0.8   -3.5      0.1    80    -2.5      0.1    0.5   2.0];
    lb = [ 0.5    0.1   0.3   -0.5      0.0    0.1   3.1       0.0    0.1   -3.9      0.0    20    -3.0      0.0    0.1   1.5];
    ub = [ 1.1    20    5.0   +0.5      0.8    8.0   3.9       0.8    15.0  -3.1      0.8    600   -2.0      0.8    13    2.5];
    fit_para = zeros(length(iv), Npixel);
    Zamide_DZ_vec = zeros(1,Npixel); % based on delta Z
    Zguan_DZ_vec = zeros(1,Npixel);
    Znoe_DZ_vec = zeros(1,Npixel);

    offs_bak_ind = find( offs <= 1 | offs >= 6 );
    % offs_bak_ind = find( abs(offs)>=5 | abs(offs)<=1 );
    offs_bak = offs(offs_bak_ind);
    fit_bak_ind = [1:4,8:13]; % Water + MT + NOE
    % fit_bak_ind = [1:4,11:13]; % Water + MT

    %% MPLF fitting
    backNum = 0;
    tic
    for i = 1:Npixel
        fprintf(repmat('\b',1,backNum));
        backNum = fprintf('Calculation process: %d/%d',i,Npixel);

        zspec = zspec_vec(:,i);

        %====================== first fitting (DS+MT+NOE) =======================
        % par_bak = lsqcurvefit(MPLF_ss,iv(fit_bak_ind),offs_bak,zspec(offs_bak_ind),lb(fit_bak_ind),ub(fit_bak_ind),OPTIONS);
        problem = createOptimProblem('lsqcurvefit','x0',iv(fit_bak_ind),'objective',R1rho_ss,...
    'lb',lb(fit_bak_ind),'ub',ub(fit_bak_ind),'xdata',offs_bak,'ydata',zspec(offs_bak_ind));

        rng default % for reproducibility
        ms = MultiStart('Display','off');
        [par_bak,~] = run(ms,problem,multistart_N);

        % % check if reach boundary
        % if ~isempty( find( par_bak==lb(fit_bak_ind) | par_bak==ub(fit_bak_ind), 1) )
        %     keyboard
        % end

        %====================== second fitting (fix par_bak) =======================
        iv_temp = iv;iv_temp(fit_bak_ind) = par_bak;
        lb_temp = lb;lb_temp(fit_bak_ind) = par_bak;
        ub_temp = ub;ub_temp(fit_bak_ind) = par_bak;
        problem = createOptimProblem('lsqcurvefit','x0',iv_temp,'objective',R1rho_ss,...
    'lb',lb_temp,'ub',ub_temp,'xdata',offs,'ydata',zspec);
        ms = MultiStart('Display','off');
        [par,~] = run(ms,problem,multistart_N);
        fit_para(:,i) = par;

        % % check if reach boundary
        % if ~isempty( find( par==lb | par==ub, 1) )
        %     keyboard
        % end

        % fitted curves
        Z_fit1 = R1rho_ss(par_bak,offs);
        Z_fit2 = R1rho_ss(par,offs);
        % Z_bak = R1rho_ss([par(1:4),par(11:13)],offs); % Water + MT
        Z_guan = R1rho_ss(par([1:4,8:16]), offs); % Water + MT + NOE + guan
        Z_amide = R1rho_ss(par([1:7,8:13]), offs); % Water + MT + NOE + amide
        Z_bak = R1rho_ss([par(1:4),par(8:13)],offs); % Water + MT + NOE

        Zamide_DZ_vec(i) = 100*max(Z_bak-Z_amide);
        Zguan_DZ_vec(i) = 100*max(Z_bak-Z_guan);
        % Znoe_DZ_vec(i) = 100*max(Z_bak-Z_noe);
    end
    fprintf(repmat('\b',1,backNum));
    fprintf('elapsed time: %.4f s\n', toc);

end