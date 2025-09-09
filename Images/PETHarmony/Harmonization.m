% multi-site SSM
clc,clear

%% load data
site1_struct = get_file_info('C:\ModelData\PET\HuaShan\NC\nii');
site2_struct = get_file_info('C:\ModelData\PET\HuaShan\MCI\nii');
site3_struct = get_file_info('C:\ModelData\PET\ADNI\NC\nii');
site4_struct = get_file_info('C:\ModelData\PET\ADNI\MCI\nii');
site1_struct = [site1_struct; site2_struct];
site2_struct = [site3_struct; site4_struct];
first_image_info = spm_vol([site1_struct(1).folder,'\',site1_struct(1).name]);
mask_name = 'C:\ModelData\PET\mask12815.nii';
mask_head=spm_vol(mask_name);
% mask_head.fname = 'mask_up0.nii';
% mask_data = spm_read_vols(spm_vol(mask_name));

[dim_x,dim_y,dim_z] = size(spm_read_vols(first_image_info));
site1_image_num = length(site1_struct);
image_num = site1_image_num+length(site2_struct);
harmonize_data = zeros(image_num,prod([dim_x,dim_y,dim_z]));

%% covariate_info
% load F:\MATLAB\MultiSite-SVD\age.mat;
% load F:\MATLAB\MultiSite-SVD\sex.mat;
% load F:\MATLAB\MultiSite-SVD\TIV.mat;
% x_age = age;
% x_tiv = TIV;
% x_age = mapminmax(age',0,1)';
% x_tiv = mapminmax(TIV',0,1)';
% x_sex = dummyvar(sex);
% load E:\SVD\cov\A.mat;
% covariate_info = A;
load C:\ModelData\huashan_adni.mat;

%% SVD
[~,subject_data,subject_score,gismatrix,mask_up0] = ssm_pca(site1_struct,site2_struct,mask_name,all_data);
% [~,subject_data,subject_score,gismatrix,mask_up0] = ssm_pca(site1_struct,site2_struct,mask_name);

%% Harmonization
Ttest = harmonize_image(subject_score,gismatrix,site1_image_num);
subject_score = subject_score(:,Ttest(:,1));
gismatrix = gismatrix(:,Ttest(:,1));
subject_score_har = harmonize_score(subject_score,site1_image_num,Ttest,1,0.05);    % Abeta 0.001  Tau 0.05

%% Write image
P = 0.001;
gismatrix_n = sum(Ttest(:,2)<P);
% site_score = subject_score(:,1:gismatrix_n);    %1:gismatrix_n);
% site_gismatrix = gismatrix(:,1:gismatrix_n);    %1:gismatrix_n);
% site_data = site_score*site_gismatrix';
% har_data = subject_data'-site_data;
% non_negative = any(har_data(:)<0);
% if ~non_negative

% % % % har_data = subject_score*gismatrix';
% % % % for site_gismatrix_i = 1:gismatrix_n
% % % %     site_data = subject_score(:,site_gismatrix_i)*gismatrix(:,site_gismatrix_i)';
% % % %     har_data = har_data-site_data;
% % % %     har_data(har_data<0.2) = har_data(har_data<0.2)+site_data(har_data<0.2);   %  A:0.2    FDG:0    Tau:1
% % % % end

har_data = subject_score_har*gismatrix'+subject_data';
% subject_score(:,2) = [];
% gismatrix(:,2) = [];
% har_data = subject_score*gismatrix';
% har_data(har_data<0) = har_data(har_data<0)+site_data(har_data<0);
harmonize_data(:,mask_up0==1) = har_data;

% for subject_i = 1:site1_image_num
%     first_image_info.fname = ['E:\IXIOASIS1\IXI' '\' site1_struct(subject_i).name];
%     spm_write_vol(first_image_info,reshape(harmonize_data(subject_i,:),dim_x,dim_y,dim_z));
% end

for subject_i = site1_image_num+1:image_num
    first_image_info.fname = ['C:\ModelData\PET\ADNI\har' '\' site2_struct(subject_i-site1_image_num).name];
    spm_write_vol(first_image_info,reshape(harmonize_data(subject_i,:),dim_x,dim_y,dim_z));
end

% first_image_info.fname = ['E:\IXIOASIS1' '\' 'MASK.nii'];
% spm_write_vol(first_image_info,reshape(mask_up0,dim_x,dim_y,dim_z));


% for site_gismatrix_i = 1:gismatrix_n
%     mask_up0(mask_up0~=0) = gismatrix(:,site_gismatrix_i)';
%     first_image_info.fname = ['E:\IXIOASIS1' '\' 'site_gismatrix_' num2str(site_gismatrix_i) '.nii'];
%     spm_write_vol(first_image_info,reshape(mask_up0,dim_x,dim_y,dim_z));
% end
%
% end