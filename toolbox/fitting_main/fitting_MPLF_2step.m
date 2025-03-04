function [Zamide_DZ_vec, Zguan_DZ_vec, fit_para] = fitting_MPLF_2step(offs, zspec_vec, multistart_N)
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
    iv = [ 0.9    0.8   2.0   0         0.05   2.0   3.5       0.1    9.0   -3.5      0.1    150   -2.5      0.01   2.0   2.0];
    lb = [ 0.5    0.4   0.5   -0.5      0.0    0.5   3.1       0.0    0.5   -3.9      0.0    30    -3.0      0.0    0.5   1.5];
    ub = [ 1.1    1.0   8.0   +0.5      0.2    15.0  3.9       0.6    20.0  -3.1      0.8    400   -2.0      0.2    20.0  2.5];
    fit_para = zeros(length(iv), Npixel);
    Zamide_DZ_vec = zeros(1,Npixel); % based on delta Z
    Zguan_DZ_vec = zeros(1,Npixel);

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

        %====================== first fitting (DS+MT) =======================
        % par_bak = lsqcurvefit(@lorentzMultipool,iv(fit_bak_ind),offs_bak,zspec(offs_bak_ind),lb(fit_bak_ind),ub(fit_bak_ind),OPTIONS);
        
        problem = createOptimProblem('lsqcurvefit','x0',iv(fit_bak_ind),'objective',@lorentzMultipool,...
    'lb',lb(fit_bak_ind),'ub',ub(fit_bak_ind),'xdata',offs_bak,'ydata',zspec(offs_bak_ind));

        rng default % for reproducibility
        ms = MultiStart('Display','off');
        [par_bak,~] = run(ms,problem,multistart_N);
        %====================== second fitting (fix par_bak) =======================
        iv_temp = iv;iv_temp(fit_bak_ind) = par_bak;
        lb_temp = lb;lb_temp(fit_bak_ind) = par_bak;
        ub_temp = ub;ub_temp(fit_bak_ind) = par_bak;
        % par = lsqcurvefit(@lorentzMultipool,iv_temp,offs,zspec,lb_temp,ub_temp,OPTIONS);
        problem = createOptimProblem('lsqcurvefit','x0',iv_temp,'objective',@lorentzMultipool,...
    'lb',lb_temp,'ub',ub_temp,'xdata',offs,'ydata',zspec);
        ms = MultiStart('Display','off');
        [par,~] = run(ms,problem,multistart_N);
        fit_para(:,i) = par;

        % % check if reach boundary
        % if ~isempty( find( par==lb | par==ub, 1) )
        %     keyboard
        % end

        % fitted curves
        % Z_bak = lorentzMultipool([par(1:4),par(11:13)],offs); % Water + MT
        % Z_guan = lorentzMultipool(par([1:4,11:16]), offs); % Water + MT + guan
        % Z_amide = lorentzMultipool(par([1:7,11:13]), offs); % Water + MT + amide
        % 
        % Zamide_DZ_vec(i) = 100*max(Z_bak-Z_amide);
        % Zguan_DZ_vec(i) = 100*max(Z_bak-Z_guan);

        Zamide_DZ_vec(i) = 100*par(5);
        Zguan_DZ_vec(i) = 100*par(14);
    end
    fprintf(repmat('\b',1,backNum));
    fprintf('elapsed time: %.4f s\n', toc);

end