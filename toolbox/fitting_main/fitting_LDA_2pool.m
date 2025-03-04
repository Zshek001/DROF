function [Zamide_DZ_vec, Zguan_DZ_vec] = fitting_LDA_2pool(offs, zspec_vec, multistart_N)

    %% LDA setting
    [nf,Npixel] = size(zspec_vec);
    fit_rgc = [-19, 0.8]; % Center range of Z-spectrum used for Lorentzian fitting, in ppm (0.5~1.5 ppm)
    fit_rge = [-19, 19]; % Edge range of Z-spectrum used for Lorentzian fitting, in ppm 

    %              1. water                    2. MT
    %      Zi   A1     G1      dw1       A2      G2      dw2   
    iv = [ 1    0.5    2       0        0.1      60      -2.5];
    lb = [ 0.5  0.02   0.3    -0.4      0.0      30      -3.0];
    ub = [ 1.1  1      100    +0.4      0.4      200     -2.0];

    % iv = iv(1:4);
    % lb = lb(1:4);
    % ub = ub(1:4);

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
        % fitted arameters
        % OPTIONS = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective','TolFun', 1e-6,'MaxIterations', 100,'Display','off');
        % par = lsqcurvefit(@lorentzMultipool,iv,offs,zspec,lb,ub,OPTIONS);

        fit_ind = get_fit_ind(zspec, offs, fit_rgc, fit_rge);
        offs_trunc = offs(fit_ind);
        zspec_trunc = zspec(fit_ind);

        problem = createOptimProblem('lsqcurvefit','x0',iv,'objective',@lorentzMultipool,...
    'lb',lb,'ub',ub,'xdata',offs_trunc,'ydata',zspec_trunc);
        
        rng default % for reproducibility
        ms = MultiStart('Display','off');
        [par,~] = run(ms,problem,multistart_N);
        % fit_para(:,i) = par;

        % % check if reach boundary
        % if ~isempty( find( par==lb | par==ub, 1) )
        %     keyboard
        % end

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