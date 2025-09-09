function file_info = spm_select_info(n,typ,mesg)
% Select the image with spm_select.m to output information about the image
%  n      - number of files [Default: Inf]
%            A single value or a range.  e.g.
%            1       - select one file
%            Inf     - select any number of files
%            [1 Inf] - select 1 to Inf files
%            [0 1]   - select 0 or 1 files
%            [10 12] - select from 10 to 12 files
%   typ    - file type [Default: 'any']
%            'any'   - all files
%            'image' - Image files (".img" and ".nii")
%                      Note that it gives the option to select individuals
%                      volumes of the images.
%            'mesh'  - Surface mesh files (".gii")
%            'xml'   - XML files
%            'mat'   - MATLAB .mat files or .txt files (assumed to contain
%                      ASCII representation of a 2D-numeric array)
%            'batch' - SPM batch files (.m or .mat)
%            'dir'   - select a directory
%            Other strings act as a filter to regexp. This means that
%            e.g. DCM*.mat files should have a typ of '^DCM.*\.mat$'
%   mesg   - a prompt [Default: 'Select files...']
%   file_info:folder name  
file_image_name = spm_select(n,typ,mesg);
file_num = size(file_image_name,1);

image_folder_struct = cell(file_num,1);
image_name_struct = cell(file_num,1);
start_char = 1;
for file_i = 1:file_num
    subject_image_name = file_image_name(file_i,:);
    image_name = strsplit(subject_image_name,',1');
    subject_char_num = length(image_name{1,1});
    image_info = strsplit(image_name{1,1},'\');
    filder_char_num = subject_char_num-length(image_info{1,end})-1;
    image_name_struct{file_i,1} = image_info{1,end};
    image_folder_struct{file_i,1} = subject_image_name(1,1:filder_char_num);
end
file_info = struct('folder',image_folder_struct,'name',image_name_struct);
end


