function [Zamide_DZ_vec, Zguan_DZ_vec, Znoe_DZ_vec, fit_para] = fitting_QDROF(offs, zspec_vec, FitParam, multistart_N)
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
    R1rho_ss = @(par,offs) lorentzMultipool_R1rho_quass(par,offs,FitParam);

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
        Z_fit = R1rho_ss(par,offs);
        Z_bak = R1rho_ss([par(1:4),par(11:13)],offs); % Water + MT
        Z_guan = R1rho_ss(par([1:4,11:16]), offs); % Water + MT + guan
        Z_amide = R1rho_ss(par([1:7,11:13]), offs); % Water + MT + amide
        Z_noe = R1rho_ss([par(1:4),par(8:13)],offs); % Water + MT + NOE

        Zamide_DZ_vec(i) = 100*max(Z_bak-Z_amide);
        Zguan_DZ_vec(i) = 100*max(Z_bak-Z_guan);
        Znoe_DZ_vec(i) = 100*max(Z_bak-Z_noe);
    end
    fprintf(repmat('\b',1,backNum));
    fprintf('elapsed time: %.4f s\n', toc);

end