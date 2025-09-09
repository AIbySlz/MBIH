function [eigenvector,eigvalue] = pca_f(matrix)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
[feature_num,subject_num] = size(matrix);
% mean_feature = mean(matrix')';
% mean_matrix = repmat(mean_feature,1,subject_num);
% matrix = matrix-mean_matrix;
[eigenvector eigvalue] = eig(matrix'*matrix);
eigvalue = diag(eigvalue);
eigvalue = eigvalue';
eigenvector = matrix*eigenvector;
gismatrix_sum_qurt = sqrt(sum(eigenvector.*eigenvector));
gismatrix_unit_matrix = repmat(gismatrix_sum_qurt,feature_num,1);
eigenvector = eigenvector./gismatrix_unit_matrix;
end

% 
% ,options
% [COEFF SCORE latent tsquared explained]=pca(featureMatrix);
% supLatent=cumsum(latent)./sum(latent);
% pa=featureMatrix*COEFF(:,1:3);
% 
% [COEFF SCORE latent tsquared explained]=pca(featureMatrix');        %矩阵太大无法分解时采用的方法
% featurevector=featureMatrix'*COEFF;
% explained=100*latent/sum(latent);
% supLatent=cumsum(latent)./sum(latent);
% supLatent(supLatent>0.8)=[];
% k=numel(supLatent);
% pa=featureMatrix*featurevector(:,1:k);