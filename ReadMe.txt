This repository contains the code for the paper titled, "JPEG-based Variable Block-Size Image Compression using CIE L*a*b* Color Space" 
published in the KSII Transactions on Internet and Information Systems, vol. 12, no. 10, pp. 5056 - 5078, 2018.
If using this repository, please cite the above paper.

################################################
HOW TO USE:-
1. Main compression and decompression code is contained in the MATLAB file 'JPEG_based_Variable_Block_Size_Compression.m'
2. Image to be compressed should be stored at the location pointed by the variable 'IPFolder' in the above MATLAB file.
3. The location can be changed by modifying the vaiable IPFolder.
4. The variable 'k10' in the main MATLAB file indicates a value of the CSF constant c which can be used as a quality setting parameter.
5. Varying the value of k10 will vary the compression achieved and the reconstructed image quality.
6. The compressed bitstream, reconstructed image and an excel file containing the CR, PSNR and other compression parameters are stored at 
   the location pointed by the variable 'OPFolder'.
7. As in the case of input folder, the output location pointed by OPFolder can be changed by modifying the variable OPFolder.





