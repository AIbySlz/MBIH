## MBIH
A toolbox for multicenter brain images harmonization  

## Requirements
* matlab 2020a or later
* SPM12 (https://www.fil.ion.ucl.ac.uk/spm/)

## Function
Eliminate site variations in multi-site roi or image data.

## HCOBE Theory
When using HCOBE to harmonize the multicenter data, the COBE was firstly used to decompose inter-center data, the decomposed heterogeneity space contained biological heterogeneity and center effect. 
To precisely remove the center effect while preserving biological heterogeneity, the center effect and biological heterogeneity should be further separated. In order to do that, we further decomposed the
 intra-center data using COBE. Since intra-center data were acquired on the same scanner and acquisition protocol, the decomposed heterogeneity space only contained biological heterogeneity. Finally, we
can accurately obtain the center effect by calculating the differences between inter-center and intra-center heterogeneity space.

## PETHarmony Method Theory
The PETHarmony method comprises a linear regression model and singular value decomposition.  Linear regression models decomposed brain images into biological factors and a subject residual profile (SRP).  
The center effect is inherent in the SRP.  To rectify the center effect, we applied singular value decomposition to further break down the SRP.  Subsequently, we harmonized the center effect in SRP by aligning 
the scores of the decomposed components to ensure uniform distribution across centers.  Finally, the harmonized brain image was derived by incorporating the biological factors into the adjusted SRP.

## Usage
   1. Rename the folder as "MBIH," subsequently include the "MBIH" folder and its subdirectories in the Matlab path, and ensure correct installation of SPM12.
   2. At the MATLAB prompt, type "MBIH" (Menu options will appear), and press "Enter".
### HCOBE Usage
   3. Click 'HCOBE'. This method is recommended for harmonizing MRI brain images.
   4. Click 'Add a New Site' and select the normalized images of a new site. If it is necessary to delete the image of last site, we can click 'Delete the Latest Site'.
   5. Select 'Reference Site'.  
   6. Select either the "Select Mask"  or the "Set Threshold"(default = 0.2 for MRI; threshold mask is determined by the routine by combining individual masks determined by evaluating the lower threshold for ach image as the fraction of the volume max).  
   7. Give a prefix of the harmonized image name (default 'h') in 'Filename Prefix'.
   8. Select the save path for harmonized images. 'Previous Dir' stands for saving to the original path; 'Custom Dir ' stands for custom saving path; 'Re-Custom Dir ' stands for re-custom saving path.
   9. Choose whether to keep the brain mask.
   10. Click 'Run' to begin harmonization of multi-site brain images.

### PETHarmony Usage
   3. Click 'PETHarmony'. This method is recommended for harmonizing PET brain images.
   4. Click 'Select Site 1 Data' and select the normalized images of a site. If you need to re-select the data, you can click "Select Site 1 Data" again.
   5. As in step 4, click "Select Site 2 Data" to add data for site 2.
   6. Set a site as the reference site in the "Reference Site" option. During the harmonization process, the data of reference site will not change, and the data of other sites will be mapped to the same distribution as the data of referencce site.
   7. Select either the "Set Threshold" or the "Select Mask". This part is to generate a mask, and the selected area of the mask will be harmonized for data.
   8. The "Covariates" option is optional. This part is used to include biological factors in the linear regression model. It is recommended to add a design matrix to the MATLAB workspace, where each column represents a factor, such as age, sex, etc. 
       Simply fill in the name of the design matrix and you are ready to add covariates.
   9. The "P Value" option is used to define the threshold for statistical analysis to determine a significant difference between two sites. The default P value is 0.05.
   10. The "Filename Prefix" option is filled with the file name prefix of the non-reference site data. The default is "h".
   11. The "Output Directories" option is used to select the storage path for harmonized images. 'Previous Dir' stands for saving to the original path; 'Custom Dir ' stands for custom saving path; 'Re-Custom Dir ' stands for re-custom saving path.
   12. The "Save Mask" option determines whether to save the mask used in the harmonization process will be saved.
   13. The "Save Site Feature" option is used to determine whether to save the components generated in the singular value decomposition will be saved.
   14. Click 'Run' to begin harmonization of multi-site brain images.

## Date
2024.11.12
