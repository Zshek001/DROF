% preprocess DICOM data
% Jianpan Huang - jphuang@hku.hk, 20230818
    % modified 20241028
clear
close all

addpath ..\toolbox\
%% load data
[filename,filepath] = uigetfile(".\invivo_data\raw*.mat",'select CEST raw data file (raw*.mat)');
load([filepath,filesep,filename],"img_all","slicenum","mat_wid","mat_hei","fov","thick","info");

%% rearrange data
sz = size(img_all);
if sz(3)>slicenum
    img4d = reshape(img_all, [sz(1),sz(2),slicenum,sz(3)/slicenum]);
end

% load freq. offs. in [ppm], include M0
if size(img4d,4)==55
    load('offs55p_Hz.mat','offs'); % user 19.txt on GE, full zspec
    offs = offs/(3*42.58);
elseif size(img4d,4)==46
    load('offs46p_Hz.mat','offs'); % user 6.txt on GE, full zspec
    offs = offs/(3*42.58);
else
    error('no compatible offset file found')
end
% [~,I]=sort(offs);figure();plot(offs(I),squeeze(img4d(154,104,I)));
% xlabel('offs [ppm]');set(gca,'XDir','reverse');title('zspec')

%% seperate cest images and m0 image
[img, offs, img_m0] = separate_z_m0_3d(img4d, offs,-300);
% -------- rule out 200ppm data point --------
idx_eff = find(abs(offs)<=130);
img = img(:,:,:,idx_eff);
offs = offs(idx_eff);
% --------------------------------------------
[xn,yn,sn,on] = size(img);

%% draw roi and segment
% manually draw roi
roi_name = 'roi_SingleSliceCEST';
temp = split(filepath,'\');
roi_filepath = strjoin(temp(1:end-1),'\');
for s = 1:sn
    img_temp = img_m0(:,:,s);
    roi_temp = draw_load_roi(roi_filepath, img_temp, [roi_name, num2str(s)], 'polygon');
    roi(:,:,s) = roi_temp; 
end

%% generate deltaB0 and Correct B0
fprintf('generating deltaB0, please hold\n')
z_corr = abs(img./(img_m0+eps).*roi); % [nx,ny,nz,nf]
db0 = generate_load_db0_3d(filepath, 'b0map', z_corr, offs, roi, 'spline');
z_corr = correct_b0_3d(z_corr, db0, offs, roi, 'spline');

%% save the preprocessed data
save([filepath,filesep,'cest_preproc.mat'],"img","img_m0",'z_corr','db0',"roi","offs");
