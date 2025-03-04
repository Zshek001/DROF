% perform fitting on invivo data
%       compare amide, guan, NOE
%       compare MTR, LDA, MPLF (MultiStart1), MPLF (MS-5), MPLF2-ss (MS-5)
% 2025-0203 for invivo data
clear
close all
warning off

addpath ..\toolbox\
addpath ..\toolbox\fitting_main\
multistart_N = 5; % times for MultiStart

OutputDir = ".\Outprocess";
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

%% three in-vivo data 
datalist = ["20241022_1","20241022_2","20241024"];
for dataidx = 1:3
    %% load data
    % load(fileName,"offs","zspec"); % m0 not included
    %     % offs: [nf,1]
    %     %   offs(1): m0
    %     % zspec: [nf, nguanf, nguank, namidef]
    % 
    % [nf,nguanf,nguank,namidef,nMTf,nNOEf] = size(zspec);
    % zspec_vec = reshape(zspec,nf,[]);
    % Npixel = size(zspec_vec,2);
    OutputDirTemp = OutputDir + datalist(dataidx) + "\";
    mkdir(OutputDirTemp);
    
    load("invivo_preproc_"+datalist(dataidx)+".mat",'z_corr','img_m0','offs','roi');
    zmap_mea = squeeze(abs(z_corr.*roi)); % [nx,ny,nf]
    [nx,ny,nf] = size(zmap_mea);
    
    roi_vec = reshape(roi,[],1);
    indnz = find(roi_vec~=0);
    zspec_vec = reshape(permute(zmap_mea,[3,1,2]),nf,[]); % [nf,Npixel]
    zspec_nz_vec = zspec_vec(:,indnz);
    
    % to be output
    amide_DZ_map = zeros(nx,ny);
    guan_DZ_map = zeros(nx,ny);
    noe_DZ_map = zeros(nx,ny);
    amide_CO_map = zeros(nx,ny);
    guan_CO_map = zeros(nx,ny);
    noe_CO_map = zeros(nx,ny);
    
    %% (1) Delta Z
    processing_DZ_methods = {@fitting_MTR,@fitting_LDA};
    output_filenames = {'DZ_MTR','DZ_LDA'};
    for i = 1:length(processing_DZ_methods)
        method = processing_DZ_methods{i};
        [Zamide_DZ_vec, Zguan_DZ_vec, Znoe_DZ_vec] = method(offs, zspec_nz_vec, multistart_N);
    
        amide_DZ_map(indnz) = Zamide_DZ_vec;
        guan_DZ_map(indnz) = Zguan_DZ_vec;
        noe_DZ_map(indnz) = Znoe_DZ_vec;
    
        save(OutputDirTemp+"output_"+output_filenames{i}+".mat",'amide_DZ_map','guan_DZ_map','noe_DZ_map','roi','offs','zmap_mea','img_m0');
    
        fprintf('saved to %s\n', "output_"+output_filenames{i}+".mat");
    end
    
    %% (2) Coefficience, MPLF
    processing_DZ_methods = {@fitting_MPLF};
    output_filenames = {'CO_MPLF'};
    
    
    method = processing_DZ_methods{1};
    [Zamide_DZ_vec, Zguan_DZ_vec, Znoe_DZ_vec,fit_paratemp] = method(offs, zspec_nz_vec, multistart_N);
    
    amide_DZ_map(indnz) = Zamide_DZ_vec;
    guan_DZ_map(indnz) = Zguan_DZ_vec;
    noe_DZ_map(indnz) = Znoe_DZ_vec;
    amide_CO_map(indnz) = fit_paratemp(5,:);
    guan_CO_map(indnz) = fit_paratemp(14,:);
    noe_CO_map(indnz) = fit_paratemp(8,:);
    
    nterm = size(fit_paratemp,1);
    fit_para = zeros(nterm,nx*ny);
    fit_para(:,indnz) = fit_paratemp;
    fit_para = reshape(fit_para,[nterm,nx,ny]);
    
    save(OutputDirTemp+"output_"+output_filenames{1}+".mat",'amide_DZ_map','guan_DZ_map','noe_DZ_map','roi','offs','fit_para',...
                                                                    'amide_CO_map','guan_CO_map','noe_CO_map');
    fprintf('saved to %s\n', "output_"+output_filenames{1}+".mat");
    
    
    %% (3) R1rho, DROF
    B0 = 3; % T
    FitParam.R1 = 1; % brain ~ 1s @ 3T
    FitParam.satpwr = 0.8; % [uT]
    FitParam.tsat = 2; % [s]
    FitParam.Magfield = 42.5764 * B0; % [Hz]
    FitParam.PeakOffset = 3.5; % used in PLOF, position to be polynomialize
    
    processing_R1rho_methods = {@fitting_DROF};
    output_filenames = {'R1rho_DROF'};
    for i = 1:length(processing_R1rho_methods)
        method = processing_R1rho_methods{i};
        [Zamide_DZ_vec, Zguan_DZ_vec, Znoe_DZ_vec, fit_paratemp] = method(offs, zspec_nz_vec, FitParam, multistart_N); 
    
        amide_DZ_map(indnz) = Zamide_DZ_vec;
        guan_DZ_map(indnz) = Zguan_DZ_vec;
        noe_DZ_map(indnz) = Znoe_DZ_vec;
        amide_CO_map(indnz) = fit_paratemp(5,:);
        guan_CO_map(indnz) = fit_paratemp(14,:);
        noe_CO_map(indnz) = fit_paratemp(8,:);
    
        nterm = size(fit_paratemp,1);
        fit_para = zeros(nterm,nx*ny);
        fit_para(:,indnz) = fit_paratemp;
        fit_para = reshape(fit_para,[nterm,nx,ny]);
    
        save(OutputDirTemp+"output_"+output_filenames{1}+".mat",'amide_DZ_map','guan_DZ_map','noe_DZ_map','roi','offs','fit_para',...
                                                                        'amide_CO_map','guan_CO_map','noe_CO_map');
        fprintf('saved to %s\n', "output_"+output_filenames{i}+".mat");
    end

end