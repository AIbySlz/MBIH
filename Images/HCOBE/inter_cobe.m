function [original_voxel_matrix, voxel_matrix,size_matrix,mask] = inter_cobe(site_file,site_image_num,mask_type,reference_site)
%   inter-site cobe
%   voxel_matrix:raw space + common space  size_matrix:image size
site_num = length(site_image_num);
image_num = sum(site_image_num);

%% Initialization
image_name = [site_file(1).folder,'\',site_file(1).name];
image_matrix = spm_read_vols(spm_vol(image_name));
size_matrix = size(image_matrix);
voxel_matrix = zeros(prod(size_matrix),image_num);
original_voxel_matrix = voxel_matrix;

%% Mask and read data
if isnumeric(mask_type)
    h_image = waitbar(0,'Image reading process (Inter-Site)');
    mask = ones(prod(size_matrix),1);
    for image_i = 1:image_num
        s = ['Image reading process (Inter-Site): ' num2str(ceil(100*image_i/image_num)) '%'];
        image_name = [site_file(image_i).folder,'\',site_file(image_i).name];
        image_matrix = spm_read_vols(spm_vol(image_name));
        voxel_matrix(:,image_i) = image_matrix(:);
        original_voxel_matrix(:,image_i) = image_matrix(:);
        subject_max = max(voxel_matrix(:,image_i));
        threshold_s = mask_type*subject_max;
        mask_s = (voxel_matrix(:,image_i) >= threshold_s) ;
        mask = mask.*mask_s ;
        waitbar(image_i/image_num,h_image,s);
    end
    close(h_image);
    disp(strcat(datestr(datetime),'-Done    ''Inter-Site: Reading Data'''));
    
else
    mask = spm_read_vols(spm_vol(mask_type));
    mask(mask~=0) = 1;
    if ~isequal(size_matrix,size(mask))
        error('The image dimensions of subjects is not equal to the dimension of the MASK');
    end
    
    mask = mask(:);
    h_image = waitbar(0,'Image reading process (Inter-Site)');
    for image_i = 1:image_num
        s = ['Image reading process (Inter-Site): ' num2str(ceil(100*image_i/image_num)) '%'];
        image_name = [site_file(image_i).folder,'\',site_file(image_i).name];
        image_matrix = spm_read_vols(spm_vol(image_name));
        voxel_matrix(:,image_i) = image_matrix(:);
        original_voxel_matrix(:,image_i) = image_matrix(:);
        mask = mask.*(voxel_matrix(:,image_i) > 0);
        waitbar(image_i/image_num,h_image,s);
    end
    close(h_image);
    disp(strcat(datestr(datetime),'-Done    ''Inter-Site: Reading Data'''));
    
end

%% Mask image
voxel_matrix(mask==0,:)=[];

%% PCA
% [~,~,~,~,explained] = pca(voxel_matrix);
% % sum_score = cumsum(explained)/100;
% % feature_num = sum(sum_score<0.7);
% explained = explained/sum(explained);
% feature_num = sum(explained>0.01);

%% COBE
Y = cell(1,site_num);
for site_i = 1:site_num
    start_image = sum(site_image_num(1:site_i-1))+1;
    end_image = sum(site_image_num(1:site_i));
    eval(['Y{' num2str(site_i) '} = voxel_matrix(:,start_image:end_image);']);
end

opts.c = [];                % number of common features
[A,B,~] = cobe(Y,opts);  % common features extraction
B = mape_b(B,reference_site);

for site_i = 1:site_num
    start_image = sum(site_image_num(1:site_i-1))+1;
    end_image = sum(site_image_num(1:site_i));
    eval(['voxel_matrix(:,start_image:end_image) = voxel_matrix(:,start_image:end_image)+A*B{1,' num2str(site_i) '};']);
end
disp(strcat(datestr(datetime),'-Done    ''Inter-Site: COBE'''));

end

