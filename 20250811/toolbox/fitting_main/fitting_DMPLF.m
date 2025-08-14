function [varargout] = fitting_DMPLF(offs, zspec_vec, multistart_N, display_on)
% INPUT
%   offs: [nf,1]
%   zspec_vec: [nf,Npixel]
% OUTPUT
%   Zamide_DZ_vec: [1,Npixel]
%   Zguan_DZ_vec: [1,Npixel]
%   fit_para: [16, Npixel]
%       Zi (1), Water (2-4), Amide (5-7), NOE (8-10), MT (11-13), Guan (14-16)

    if ~exist('display_on','var') | isempty(display_on)
        display_on = 1;
    end
    rng default % for reproducibility

    %% MPLF setting
    [nf,Npixel] = size(zspec_vec);

    %             1. Water              2. Amide               3. NOE                 4. MT                  5. Guan
    %      Zi     A1    G1    dw1       A2     G2    dw2       A3     G3    dw3       A4     G4    dw4       A5     G5    dw5
    iv = [ 0.9    0.9   2.3   0         0.05   2.0   3.5       0.1    9.0   -3.5      0.1   150    -2.5      0.01   2.0   2.0];
    lb = [ 0.5    0.1   0.5   -0.5      0.0    0.5   3.1       0.0    0.5   -3.9      0.0    30    -3.0      0.0    0.5   1.5];
    ub = [ 1.1    1     8.0   +0.5      0.25   15.0  3.9       0.6   20.0   -3.1      0.8   400    -2.0      0.25   30    2.5];
    
    fit_para = zeros(length(iv), Npixel);

    offs_bak_ind = find( offs <= 1 | offs >= 8 );
    offs_bak = offs(offs_bak_ind);
    fit_bak_ind = [1:4,8:13]; % Water + MT + NOE

    %% MPLF fitting
    tic
    runMS = @(prob, n) run(MultiStart('Display','off'), prob, n);
    backNum = 0;
    fprintf('\t processing DMPLF: ');
    for i = 1:Npixel
        if display_on == 1
            fprintf(repmat('\b',1,backNum));
            backNum = fprintf('%d/%d',i,Npixel);
        end
        zspec = zspec_vec(:,i);

        %====================== first fitting (DS+MT) =======================
        problem = createOptimProblem('lsqcurvefit','x0',iv(fit_bak_ind),'objective',@lorentzMultipool,...
    'lb',lb(fit_bak_ind),'ub',ub(fit_bak_ind),'xdata',offs_bak,'ydata',zspec(offs_bak_ind));

        [par_bak, ~] = runMS(problem, multistart_N);
        %====================== second fitting (fix par_bak) =======================
        iv_temp = iv;iv_temp(fit_bak_ind) = par_bak;
        lb_temp = lb;lb_temp(fit_bak_ind) = par_bak;
        ub_temp = ub;ub_temp(fit_bak_ind) = par_bak;
        problem = createOptimProblem('lsqcurvefit','x0',iv_temp,'objective',@lorentzMultipool,...
    'lb',lb_temp,'ub',ub_temp,'xdata',offs,'ydata',zspec);

        [par, ~] = runMS(problem, multistart_N);
        fit_para(:,i) = par;
    end
    if display_on == 1
        fprintf(repmat('\b',1,backNum));
        fprintf('elapsed time: %.4f s\n', toc);
    end

    switch nargout
        case 1
            varargout{1} = fit_para;
        case 3
            Zamide_DZ_vec = 100*fit_para(5,:);
            Zguan_DZ_vec = 100*fit_para(14,:);
            Znoe_DZ_vec = 100*fit_para(8,:);
            varargout = {Zamide_DZ_vec, Zguan_DZ_vec, Znoe_DZ_vec};
        otherwise
            error('Unsupported output parameters number');
    end

end