function [Zamide_DZ_vec, Zguan_DZ_vec, Znoe_DZ_vec] = fitting_MTR_par(offs, zspec_vec, ~, displayflag)

 
    if ~exist('displayflag','var') || isempty(displayflag)
        displayflag = 1;
    end
    %% MTR setting
    [nf,Npixel] = size(zspec_vec);

    [~, ind_amide] = min(abs(offs - 3.5));
    [~, ind_amide_n] = min(abs(offs + 3.5));
    [~, ind_guan] = min(abs(offs - 2));
    [~, ind_guan_n] = min(abs(offs + 2));

    Zamide_DZ_vec = zeros(1,Npixel);
    Zguan_DZ_vec = zeros(1,Npixel);
    Znoe_DZ_vec = zeros(1,Npixel);

    %% MTR analysis
    tic
    if displayflag == 1
        fprintf('\t processing MTRasym, ');
    end
    for i = 1:Npixel
       
        zspec = zspec_vec(:,i);

        Zamide_DZ_vec(i) = 100*(zspec(ind_amide_n) - zspec(ind_amide));
        Zguan_DZ_vec(i) = 100*(zspec(ind_guan_n) - zspec(ind_guan));
        Znoe_DZ_vec(i) = 100*(zspec(ind_amide) - zspec(ind_amide_n));
    end
    
    if displayflag == 1
        fprintf('elapsed time: %.4f s\n', toc);
    end
end