function [covariate_data,subject_score,gismatrix] = ssm_pca_radiomics(site1_struct,site2_struct, covariate_info)
% Data input from two centers,and corrects the multi-centric effect of two centers
% Data information
% site1_struct，site2_struct：The image structure with two centers
% mask_type： Mask Image name or threshold

%% Read excel data
[site1_data,~,~] = xlsread(site1_struct);
[site2_data,~,~] = xlsread(site2_struct);
subject_data = [site1_data site2_data];
[~, site1_image_num] = size(site1_data);
[~, image_num] = size(subject_data);

%% GLM
covariate_data = zeros(size(subject_data));
if nargin == 3
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

