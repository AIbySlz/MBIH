function voxel_matrix = intra_cobe(site_file,site_image_num,mask,site_i)
%   intra-site cobe
%   voxel_matrix:common space
center1_Num = ceil(site_image_num/2);

%% Initialization
voxel_matrix = zeros(sum(mask~=0),site_image_num);

%% Read data
s = ['Image reading process (Intra-Site' num2str(site_i) ')'];
h_image = waitbar(0,s);
for image_i = 1:site_image_num
    s = ['Image reading process (Intra-Site' num2str(site_i) '):' num2str(ceil(100*image_i/site_image_num)) '%'];
    image_name = [site_file(image_i).folder,'\',site_file(image_i).name];
    image_matrix = spm_read_vols(spm_vol(image_name));
    image_matrix = image_matrix(:);
    voxel_matrix(:,image_i) = image_matrix(mask~=0,:);
    waitbar(image_i/site_image_num,h_image,s);
end
close(h_image);
s = ['-Done    ''Intra-Site' num2str(site_i) ': Reading Data'''];
disp(strcat(datestr(datetime),s));

%% PCA
% [~,~,~,~,explained] = pca(voxel_matrix);
% % sum_score = cumsum(explained)/100;
% % feature_num = sum(sum_score<0.7);
% explained = explained/sum(explained);
% feature_num = sum(explained>0.01);

%% COBE
Y = cell(1,2);
Y{1} = voxel_matrix(:,1:center1_Num);
Y{2} = voxel_matrix(:,center1_Num+1:site_image_num);

opts.c = [];
[A,B,~]=cobe(Y,opts);  % common features extraction

B = mape_b(B,1);
voxel_matrix(:,1:center1_Num) = A*B{1,1};
voxel_matrix(:,center1_Num+1:site_image_num) = A*B{1,2};

s = ['-Done    ''Intra-Site' num2str(site_i) ': COBE'''];
disp(strcat(datestr(datetime),s));

end

