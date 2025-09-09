function [original_data,covariate_data,subject_score,gismatrix,mask] = ssm_pca(site1_struct,site2_struct,mask_type,covariate_info)
% Data input from two centers,and corrects the multi-centric effect of two centers
% Data information
% site1_struct，site2_struct：The image structure with two centers
% mask_type： Mask Image name or threshold

%% Read image information
site1_first_name = [site1_struct(1).folder,'\',site1_struct(1).name];
site2_first_name = [site2_struct(1).folder,'\',site2_struct(1).name];
if ~isequal(size(spm_read_vols(spm_vol(site1_first_name))),size(spm_read_vols(spm_vol(site2_first_name))))
    error('The image dimensions of the two centers are not equal');
end
[dim_x,dim_y,dim_z] = size(spm_read_vols(spm_vol(site1_first_name)));
voxel_num = dim_x*dim_y*dim_z;
site1_image_num = length(site1_struct);
image_num = site1_image_num+length(site2_struct);

%% Initialization
subject_data = zeros(voxel_num,image_num);
mask = ones(voxel_num,1);

%% Read images
h_image = waitbar(0,'Image reading process');
for subject_i = 1:site1_image_num
    s = ['Image reading process: ' num2str(ceil(100*subject_i/image_num)) '%'];
    subject_name = [site1_struct(subject_i).folder,'\',site1_struct(subject_i).name];
    subject_voxel = spm_read_vols(spm_vol(subject_name));
    subject_data(:,subject_i) = subject_voxel(:);
    waitbar(subject_i/image_num,h_image,s);
end

for subject_i = site1_image_num+1:image_num
    s = ['Image reading process: ' num2str(ceil(100*subject_i/image_num)) '%'];
    subject_name = [site2_struct(subject_i-site1_image_num).folder,'\',site2_struct(subject_i-site1_image_num).name];
    subject_voxel = spm_read_vols(spm_vol(subject_name));
    subject_data(:,subject_i) = subject_voxel(:);
    waitbar(subject_i/image_num,h_image,s);
end

close(h_image);
disp(strcat(datestr(datetime),'-Done    ''Reading Data'''));
%% Read in or generate image mask
if isnumeric(mask_type)
    
    for subject_i = 1:image_num
        subject_max = max(subject_data(:,subject_i));
        threshold_s = mask_type*subject_max;
        mask_s = (subject_data(:,subject_i) >= threshold_s) ;
        mask = mask | mask_s;
        % mask = mask.*mask_s ;
    end
    
else
    
    mask = spm_read_vols(spm_vol(mask_type));
    mask(mask~=0) = 1;
    % if ~isequal([dim_x,dim_y,dim_z],size(mask))
    %     error('The image dimensions of subjects is not equal to the dimension of the MASK');
    % end
    % 
    % mask = mask(:);
    % for subject_i = 1:image_num
    %     mask = mask.*(subject_data(:,subject_i) > 0);
    % end
    
end

%% Mask images
original_data = subject_data';     % 新添
subject_data(mask==0,:)=[];

%% GLM
covariate_data = zeros(size(subject_data));
if nargin == 4
    h_image = waitbar(0,'Generalized linear model process: 0%');
    voxel_num = size(subject_data,1);
    x_intercept = ones(image_num,1);
    x_site1 = zeros(image_num,1);
    x_site1(1:site1_image_num,1) = 1;
    x_site2 = zeros(image_num,1);
    x_site2(site1_image_num+1:image_num,1) = 1;
    covariate_num = size(covariate_info,2);
    X = [covariate_info x_site1 x_site2 x_intercept];
    
    process_start = 0;
    for voxel_i = 1:voxel_num
        Y = subject_data(voxel_i,:);
        beta = pinv(X'*X)*X'*Y';
        covariate_data(voxel_i,:) = beta(1:covariate_num,1)'*X(:,1:covariate_num)';
        process_ing = ceil(100*voxel_i/voxel_num);
        if process_ing>process_start
            s = ['Generalized linear model process: ' num2str(process_ing) '%'];
            waitbar(process_ing/100,h_image,s);
            process_start = process_ing;
        end
    end
    subject_data = subject_data-covariate_data;
    close(h_image);
end

%% Compute sub-by-sub Covariance matrix
h_image = waitbar(0,'Singular value decomposition process');
subject_covar = subject_data'*subject_data;        %devs'*devs
s=['Singular value decomposition process: 20%'];
waitbar(0.2,h_image,s);

%% Compute Eigenvector (SSF) and Eigenvalue of Covariance matrix
[eigenmatrix eigenvalue] = eig(subject_covar);
s=['Singular value decomposition process: 40%'];
waitbar(0.4,h_image,s);
eigenvalue = diag(eigenvalue);
eigenvalue(eigenvalue<0)=0;
sqrt_eigenval = sqrt(eigenvalue');
sqrt_eigenvalmat = repmat(sqrt_eigenval,image_num,1);
subject_score = eigenmatrix.*sqrt_eigenvalmat;
s=['Singular value decomposition process: 60%'];
waitbar(0.6,h_image,s);
gismatrix = subject_data*eigenmatrix;
s=['Singular value decomposition process: 80%'];
waitbar(0.8,h_image,s);
gismatrix_sum_qurt = sqrt(sum(gismatrix.*gismatrix));
gismatrix_unit_matrix = repmat(gismatrix_sum_qurt,size(gismatrix,1),1);
gismatrix = gismatrix./gismatrix_unit_matrix;
s=['Singular value decomposition process: 100%'];
waitbar(1,h_image,s);
close(h_image);
disp(strcat(datestr(datetime),'-Done    ''Singular Value Decomposition'''));
end

