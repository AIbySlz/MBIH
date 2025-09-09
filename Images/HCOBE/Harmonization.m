%The COBE method was used to coordinate the data between the three centers
clc
clear

%% load data
file_1 = get_file_info('');      %enter the folder where center1 resides
file_2 = get_file_info('');      %enter the folder where center2 resides
file_3 = get_file_info('');      %enter the folder where center3 resides
mask = 'D:\IHEPComputer\D\MultiSiteProgramHSVD\brain_mask2mm.nii';
site_file = [file_1;file_2;file_3];
site_image_num =[length(file_1);length(file_2);length(file_3)];
[original_voxel_matrix, final_matrix,size_matrix,mask] = hcobe(site_file,site_image_num,mask,3);

image_num = sum(site_image_num);
harmonize_data = zeros(size(mask));
h_image = waitbar(0,'Writing inprocess');
for subject_i = 1:site_image_num(1)
    s=['Writing inprocess:' num2str(ceil(100*subject_i/image_num)) '%'];
    image_name = [site_file(subject_i).folder,'\',site_file(subject_i).name];
    image_header = spm_vol(image_name);
    image_header.fname = ['D:\IHEPComputer' '\h' site_file(subject_i).name];
    harmonize_data = original_voxel_matrix(:,subject_i);
    harmonize_data(mask~=0,:) = final_matrix(:,subject_i);
    spm_write_vol(image_header,reshape(harmonize_data,size_matrix));
    waitbar(subject_i/image_num,h_image,s);
end

for subject_i = site_image_num(1)+1:site_image_num(1)+site_image_num(2)
    s=['Writing inprocess:' num2str(ceil(100*subject_i/image_num)) '%'];
    image_name = [site_file(subject_i).folder,'\',site_file(subject_i).name];
    image_header = spm_vol(image_name);
    image_header.fname = ['E:\MATLAB\NEW\ADNIPhilipsAchievaAD3T' '\h' site_file(subject_i).name];
    harmonize_data(mask~=0,:) = final_matrix(:,subject_i);
    spm_write_vol(image_header,reshape(harmonize_data,size_matrix));
    waitbar(subject_i/image_num,h_image,s);
end

for subject_i = site_image_num(1)+site_image_num(2)+1:image_num
    s=['Writing inprocess:' num2str(ceil(100*subject_i/image_num)) '%'];
    image_name = [site_file(subject_i).folder,'\',site_file(subject_i).name];
    image_header = spm_vol(image_name);
    image_header.fname = ['E:\MATLAB\NEW\OASIS3TrioTimAD3Tselect' '\h' site_file(subject_i).name];
    harmonize_data(mask~=0,:) = final_matrix(:,subject_i);
    spm_write_vol(image_header,reshape(harmonize_data,size_matrix));
    waitbar(subject_i/image_num,h_image,s);
end
close(h_image);
