function subject_score = harmonize_score(subject_score,site1_image_num,Ttest,reference_site,P)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
gismatrix_n = sum(Ttest(:,2)<=P);
image_num = size(subject_score,1);

switch reference_site
    case 1
        for site_gismatrix_i = 1:gismatrix_n
            site1_score = subject_score(1:site1_image_num,site_gismatrix_i);
            site2_score = subject_score(site1_image_num+1:image_num,site_gismatrix_i);
            site1_mean = mean(site1_score);
            site1_std = std(site1_score);
            site2_mean = mean(site2_score);
            site2_std = std(site2_score);
            site2_score_standard = (site2_score-site2_mean)/site2_std;
            site2_score = site2_score_standard*site1_std+site1_mean;
            subject_score(site1_image_num+1:image_num,site_gismatrix_i) = site2_score;
        end
        
    case 2
        for site_gismatrix_i = 1:gismatrix_n
            site1_score = subject_score(1:site1_image_num,site_gismatrix_i);
            site2_score = subject_score(site1_image_num+1:image_num,site_gismatrix_i);
            site1_mean = mean(site1_score);
            site1_std = std(site1_score);
            site2_mean = mean(site2_score);
            site2_std = std(site2_score);
            site1_score_standard = (site1_score-site1_mean)/site1_std;
            site1_score = site1_score_standard*site2_std+site2_mean;
            subject_score(1:site1_image_num,site_gismatrix_i) = site1_score;
        end
        
end
disp(strcat(datestr(datetime),'-Done    ''Site Subject Score Mapping'''));
end

