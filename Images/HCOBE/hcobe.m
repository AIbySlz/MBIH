function [original_voxel_matrix, final_matrix,size_matrix,mask] = hcobe(site_file,site_image_num,mask_type,reference_site)
% HCOBE
% final_matrix：corrected image size_matrix：image size 
site_num = length(site_image_num);
image_num = sum(site_image_num);

%% Inter-site cobe
[original_voxel_matrix, final_matrix,size_matrix,mask] = inter_cobe(site_file,site_image_num,mask_type,reference_site);

%% Intra-site cobe and hcobe
for site_i = 1:site_num
    start_image = sum(site_image_num(1:site_i-1))+1;
    end_image = sum(site_image_num(1:site_i));
    single_site_file = site_file(start_image:end_image);
    single_site_image_num = site_image_num(site_i);
    site_matrix = intra_cobe(single_site_file,single_site_image_num,mask,site_i);
    final_matrix(:,start_image:end_image) = final_matrix(:,start_image:end_image)-site_matrix;
end

end

