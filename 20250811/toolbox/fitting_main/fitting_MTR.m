function [Zamide_DZ_vec, Zguan_DZ_vec, Znoe_DZ_vec] = fitting_MTR(offs, zspec_vec, ~)

 
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
    fprintf('\t processing MTRasym, ');
    for i = 1:Npixel
       
        zspec = zspec_vec(:,i);

        Zamide_DZ_vec(i) = 100*(zspec(ind_amide_n) - zspec(ind_amide));
        Zguan_DZ_vec(i) = 100*(zspec(ind_guan_n) - zspec(ind_guan));
        Znoe_DZ_vec(i) = 100*(zspec(ind_amide) - zspec(ind_amide_n));
    end
    
    fprintf('elapsed time: %.4f s\n', toc);
end