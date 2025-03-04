function [Zamide_DZ_vec_bes, Zguan_DZ_vec_bes, Zamide_DZ_vec_wor, Zguan_DZ_vec_wor, fit_para] = fitting_MPLF_beworst(offs, zspec_vec, multistart_N)
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
    Zamide_DZ_vec_bes = zeros(1,Npixel); % based on delta Z
    Zguan_DZ_vec_bes = zeros(1,Npixel);
    Zamide_DZ_vec_wor = zeros(1,Npixel);
    Zguan_DZ_vec_wor = zeros(1,Npixel);

    %% MPLF fitting
    backNum = 0;
    tic
    for i = 1:Npixel
        fprintf(repmat('\b',1,backNum));
        backNum = fprintf('Calculation process: %d/%d',i,Npixel);

        zspec = zspec_vec(:,i);

        problem = createOptimProblem('lsqcurvefit','x0',iv,'objective',@lorentzMultipool,...
    'lb',lb,'ub',ub,'xdata',offs,'ydata',zspec);

        rng default % for reproducibility
        ms = MultiStart('Display','off','OutputFcn', @record_max_fval_params);
        [par,~] = run(ms,problem,multistart_N);
        Zamide_DZ_vec_bes(i) = 100*par(5);
        Zguan_DZ_vec_bes(i) = 100*par(14);

        par = record_max_fval_params();
        Zamide_DZ_vec_wor(i) = 100*par(5);
        Zguan_DZ_vec_wor(i) = 100*par(14);
    end
    fprintf(repmat('\b',1,backNum));
    fprintf('elapsed time: %.4f s\n', toc);

end