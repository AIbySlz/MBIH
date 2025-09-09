function Ttest = harmonize_image(subject_score,gismatrix,site1_image_num)
% The subject scores were statistically analyzed by T test
% subject_score:N×N    gismatrix:M×N     site1_image_num: Number of images of center 1
%% Read data
subject_num = size(gismatrix,2);

%% Initialization
Ttest = ones(subject_num,2);
Ttest(:,1) = (1:subject_num)';

%% T_test
for gis_i = 1:subject_num
    site1_socre = subject_score(1:site1_image_num,gis_i);
    site2_socre = subject_score(site1_image_num+1:subject_num,gis_i);
    [~,Ttest(gis_i,2)] = ttest2(site1_socre,site2_socre);
end

%% The GIS is in descending order according to P value
Ttest = sortrows(Ttest,2);

end

