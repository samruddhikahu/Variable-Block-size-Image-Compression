% This program contains the main code for the compression and reconstruction of an uncompressed image. 
% All the other functions in this repository are used in this code.
%Input: - An uncompressed image in .bmp or .tif file format stored at the
%         location represented by the variable IPFolder in this code.
%         Quality factor in the form of CSF constant c represented by the variable k10 in this code.
%         Changing the above two parameters will achieve different compression
%         ratios and reconstructed image qualities.
%
%Output:- A reconstructed image file in .bmp file format and 
%         the compressed bitstream stored in a .txt file stored at the
%         location represented by the variable OPFolder in this code.
%         Compression ratio, MSE, PSNR, SSIM, Compressed File size in kBs, 
%         time required for compression and decompression of the input image 
%         are stored in an excel file 'Data.xls' in the OPFolder.

clc;
clear all;
close all;

tic;
IPFolder = fullfile('D:\MATLAB','Input'); % the image(s) to be compressed should be stored at this path.
OPFolder = fullfile('D:\MATLAB','Output'); %the reconstructed image is stored at this path.
if ~exist(OPFolder, 'dir')
    mkdir(OPFolder)
end
OPXFileName = fullfile(OPFolder,'Data.xls'); % the compression parameters are stored in this excel sheet which can be found at the Output folder.

const = -2:0.01:1;
ind = find(const==0);
const(ind) = [];
for k10 = 1%:length(const) % the value of k10 decides the compression and the reconstructed quality achieved.
    %k10 = 56;
    Fname = quantization_matrix_adaptive_ksii(const(k10));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% JPEG CODING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for k10 = 1:11
    Tstrt = tic;
    clear A Alab Axyz CY C1 CCb CCr B I Ig Cres Crgb CYrec CCbrec CCbrs CCrrec CCrrs;
    clear var_mat mn_mat var_mat1 mn_mat1 Hist Idx lbl_y lbl_arr_y lbl_b lbl_arr_b lbl_r lbl_arr_r;
    clear Bout_acy Bout_dcy Seq_acy Seq_acyr Seq_dcy Seq_dcyr Code_arr_y dc_y dc_yr;
    clear Bout_acb Bout_dcb Seq_acb Seq_acbr Seq_dcb Seq_dcbr Code_arr_b dc_b dc_br;
    clear Bout_acr Bout_dcr Seq_acr Seq_acrr Seq_dcr Seq_dcrr Code_arr_r dc_r dc_rr;
    clear Hist Idx;
    clc;
    close all;
    
    IPBaseFileName = sprintf('Lena.bmp'); % the name of the image to be compressed should be entered here.
    IPFullFileName = fullfile(IPFolder,IPBaseFileName);
    A = imread(IPFullFileName);
    imshow(A);
    title('Original Image');
    B = rgb2gray(A);
    %[Axyz, Alab] = rgbtoxyzLab(A);
    Alab = rgb2lab(A);
    %Alab = rgb2ycbcr(A);
    CY = Alab(:,:,1);
    CCb = Alab(:,:,2);
    CCr = Alab(:,:,3);
    
    % load('Nmat.mat');
    % load('Cnmat.mat');
    % load('DC.mat');
    % load('Table.mat');
    
    [m, n] = size(B);
    sz = [m, n];
    toc
    
    %%%%%%%%%%%%%%%%% CODING THE Y COMPONENT %%%%%%%%%%%%%%%%%%%%%%%
    tic;
    mn_mat = [];
    var_mat = [];
    for i = 1:8:m-7
        mn_mat1 = [];
        var_mat1 = [];
        for j = 1:8:n-7
            b = CY(i:i+7,j:j+7);
            b = double(b);
            % Level shifting by -128
            bshifted = b-50;
            mn = mean(mean(bshifted));
            var = mean(mean((bshifted-mn).^2));
            mn_mat1 = [mn_mat1 mn];
            var_mat1 = [var_mat1 var];
        end
        mn_mat = [mn_mat;
            mn_mat1];
        var_mat = [var_mat;
            var_mat1];
    end
    %th = std(std(var_mat)); %/2;
    %Tvar = mean(mean(var_mat))/2;
    %Tmn = abs(mean(mean(mn_mat)));
    Tmn = std(std(mn_mat));
    Tvar =std(std(var_mat));
    toc
    
    % Forming a histogram based on mean and variance of each block:
    tic;
    s = 0;
    cnt = 1;
    %Idx = zeros(m/8,n/8);
    for i = 1:(m/8)
        for j = 1:(n/8)
            flag = 'n';
            if(s>0)
                [l l1]= size(Hist);
                for k = l:-1:1
                    dv = abs(var_mat(i,j) - Hist(k,3));
                    dm = abs(mn_mat(i,j) - Hist(k,2));
                    if((dv<Tvar)&&(dm<Tmn))
                        Hist(k,4) = Hist(k,4) + 1;
                        Idx(i,j) = k;
                        flag = 'm';
                        break;
                    end
                end
            end
            if(flag=='n')
                s = s + 1;
                Hist(s,1) = s;
                Hist(s,2) = mn_mat(i,j);
                Hist(s,3) = var_mat(i,j);
                Hist(s,4) = 1;
                Idx(i,j) = s;
            end
        end
    end
    toc
    
    % Grouping and coding of the image blocks:
    tic;
    [Seq_acy, dc_y, lbl_y, lbl_arr_y] = group_and_code(CY, Idx, Fname, 1);
    toc
    
    % Encoding the differential DC coefficients using Arithmetic Coding:
    tic;
    % DC Coefficient Sequence Generation:
    tic;
    Seq_dcy = dc_seq(dc_y);
    toc
    
    % Encoding using Arithmetic QM coding:
    tic;
    Bout_dcy = QMcoder(Seq_dcy);
    toc
    toc
    
    % Encoding AC Sequence using Arithmetic QM coding:
    tic;
    Bout_acy = QMcoder(Seq_acy);
    toc
    
    % Encoding lbl_arr_y using Exponential Golomb Coding:
    tic;
    Code_arr_y = [];
    for i = 1:length(lbl_arr_y)
        clear bits;
        bits = enc_golomb(lbl_arr_y(i),0);
        Code_arr_y = [Code_arr_y bits];
    end
    toc
    
    
    %%%%%%%%%%%%%%%%%%%%%% Coding the Cb Component %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Subsampling ratio used 4:2:0.
    tic;
    clear CCbs;
    CCbs = zeros(m/2,n/2);
    k = 1;
    for i = 1:2:m
        l = 1;
        for j = 1:2:n
            G = CCb(i:i+1,j:j+1);
            Avg = (sum(sum(G)))/4;
            CCbs(k,l) = Avg;
            l = l + 1;
        end
        k = k + 1;
    end
    toc
    
    tic;
    m1 = m/2;
    n1 = n/2;
    clear mn_mat var_mat mn_mat1 var_mat1;
    mn_mat = [];
    var_mat = [];
    for i = 1:8:m1-7
        mn_mat1 = [];
        var_mat1 = [];
        for j = 1:8:n1-7
            b = CCbs(i:i+7,j:j+7);
            b = double(b);
            % Level shifting by -128
            bshifted = b;
            mn = mean(mean(bshifted));
            var = mean(mean((bshifted-mn).^2));
            mn_mat1 = [mn_mat1 mn];
            var_mat1 = [var_mat1 var];
        end
        mn_mat = [mn_mat;
            mn_mat1];
        var_mat = [var_mat;
            var_mat1];
    end
    %th = std(std(var_mat)); %/2;
    %Tvar = mean(mean(var_mat))/2;
    %Tmn = abs(mean(mean(mn_mat)));
    Tmn = std(std(mn_mat));
    Tvar =std(std(var_mat));
    toc
    
    % Forming a histogram based on mean and variance of each block:
    tic;
    s = 0;
    cnt = 1;
    clear Hist Idx;
    %Idx = zeros(m/8,n/8);
    for i = 1:(m1/8)
        for j = 1:(n1/8)
            flag = 'n';
            if(s>0)
                [l, l1]= size(Hist);
                for k = l:-1:1
                    dv = abs(var_mat(i,j) - Hist(k,3));
                    dm = abs(mn_mat(i,j) - Hist(k,2));
                    if((dv<Tvar)&&(dm<Tmn))
                        Hist(k,4) = Hist(k,4) + 1;
                        Idx(i,j) = k;
                        flag = 'm';
                        break;
                    end
                end
            end
            if(flag=='n')
                s = s + 1;
                Hist(s,1) = s;
                Hist(s,2) = mn_mat(i,j);
                Hist(s,3) = var_mat(i,j);
                Hist(s,4) = 1;
                Idx(i,j) = s;
            end
        end
    end
    toc
    
    % Grouping and coding of the image blocks:
    tic;
    [Seq_acb, dc_b, lbl_b, lbl_arr_b] = group_and_code(CCbs, Idx, Fname, 2);
    toc
    
    % Encoding the differential DC coefficients using Arithmetic Coding:
    tic;
    % DC Coefficient Sequence Generation:
    tic;
    Seq_dcb = dc_seq(dc_b);
    toc
    
    % Encoding using Arithmetic QM coding:
    tic;
    Bout_dcb = QMcoder(Seq_dcb);
    toc
    toc
    
    % Encoding AC Sequence using Arithmetic QM coding:
    tic;
    Bout_acb = QMcoder(Seq_acb);
    toc
    
    % Encoding lbl_arr_y using Exponential Golomb Coding:
    tic;
    Code_arr_b = [];
    for i = 1:length(lbl_arr_b)
        clear bits;
        bits = enc_golomb(lbl_arr_b(i),0);
        Code_arr_b = [Code_arr_b bits];
    end
    toc
    
    %%%%%%%%%%%%%%%%%%%%%% Coding the Cr Component %%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Subsampling ratio used 4:2:0.
    clear CCrs;
    CCrs = zeros(m/2,n/2);
    k = 1;
    for i = 1:2:m
        l = 1;
        for j = 1:2:n
            G = CCr(i:i+1,j:j+1);
            Avg = (sum(sum(G)))/4;
            CCrs(k,l) = Avg;
            l = l + 1;
        end
        k = k + 1;
    end
    
    tic;
    m1 = m/2;
    n1 = n/2;
    clear mn_mat var_mat mn_mat1 var_mat1;
    mn_mat = [];
    var_mat = [];
    for i = 1:8:m1-7
        mn_mat1 = [];
        var_mat1 = [];
        for j = 1:8:n1-7
            b = CCrs(i:i+7,j:j+7);
            b = double(b);
            % Level shifting by -128
            bshifted = b;
            mn = mean(mean(bshifted));
            var = mean(mean((bshifted-mn).^2));
            mn_mat1 = [mn_mat1 mn];
            var_mat1 = [var_mat1 var];
        end
        mn_mat = [mn_mat;
            mn_mat1];
        var_mat = [var_mat;
            var_mat1];
    end
    %th = std(std(var_mat)); %/2;
    %Tvar = mean(mean(var_mat))/2;
    %Tmn = abs(mean(mean(mn_mat)));
    Tmn = std(std(mn_mat));
    Tvar =std(std(var_mat));
    toc
    
    % Forming a histogram based on mean and variance of each block:
    tic;
    s = 0;
    cnt = 1;
    clear Hist Idx;
    %Idx = zeros(m/8,n/8);
    for i = 1:(m1/8)
        for j = 1:(n1/8)
            flag = 'n';
            if(s>0)
                [l l1]= size(Hist);
                for k = l:-1:1
                    dv = abs(var_mat(i,j) - Hist(k,3));
                    dm = abs(mn_mat(i,j) - Hist(k,2));
                    if((dv<Tvar)&&(dm<Tmn))
                        Hist(k,4) = Hist(k,4) + 1;
                        Idx(i,j) = k;
                        flag = 'm';
                        break;
                    end
                end
            end
            if(flag=='n')
                s = s + 1;
                Hist(s,1) = s;
                Hist(s,2) = mn_mat(i,j);
                Hist(s,3) = var_mat(i,j);
                Hist(s,4) = 1;
                Idx(i,j) = s;
            end
        end
    end
    toc
    
    % Grouping and coding of the image blocks:
    tic;
    [Seq_acr, dc_r, lbl_r, lbl_arr_r] = group_and_code(CCrs, Idx, Fname, 2);
    toc
    
    % Encoding the differential DC coefficients using Arithmetic Coding:
    tic;
    % DC Coefficient Sequence Generation:
    tic;
    Seq_dcr = dc_seq(dc_r);
    toc
    
    % Encoding using Arithmetic QM coding:
    tic;
    Bout_dcr = QMcoder(Seq_dcr);
    toc
    toc
    
    % Encoding AC Sequence using Arithmetic QM coding:
    tic;
    Bout_acr = QMcoder(Seq_acr);
    toc
    
    % Encoding lbl_arr_y using Exponential Golomb Coding:
    tic;
    Code_arr_r = [];
    for i = 1:length(lbl_arr_r)
        clear bits;
        bits = enc_golomb(lbl_arr_r(i),0);
        Code_arr_r = [Code_arr_r bits];
    end
    toc
    
    Tcomp(k10) = toc(Tstrt)
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% JPEG DECODING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%% Decoding the Y Component %%%%%%%%%%%%
    
    Tstrt1 = tic;
    
    %%% Decoding DC:
    % Decoding using Arithmetic QM coding;
    tic;
    Seq_dcyr = QMdecoder(Bout_dcy, length(Seq_dcy));
    toc
    
    % Reconstructing DC coefficients using Seq_dcr:
    tic;
    dc_yr = dc_seq_re1(Seq_dcyr, length(dc_y));
    toc
    
    %%% Decoding AC:
    % Decoding using Arithmetic QM coding;
    tic;
    Seq_acyr = QMdecoder(Bout_acy, length(Seq_acy));
    toc
    
    % Decoding the variable sized blocks and placing them:
    tic;
    %CYrec = zeros(size(CY));
    CYrec = degroup_and_decode(Seq_acyr, dc_yr, lbl_y, Fname, n, 1);
    toc
    
    %%%%%%%%%%%%%%%%%%% Decoding the Cb Component %%%%%%%%%%%%%%%%%%%%%%
    
    %%% Decoding DC:
    % Decoding using Arithmetic QM coding;
    tic;
    Seq_dcbr = QMdecoder(Bout_dcb, length(Seq_dcb));
    toc
    
    % Reconstructing DC coefficients using Seq_dcr:
    tic;
    dc_br = dc_seq_re1(Seq_dcbr, length(dc_b));
    toc
    
    %%% Decoding AC:
    % Decoding using Arithmetic QM coding;
    tic;
    Seq_acbr = QMdecoder(Bout_acb, length(Seq_acb));
    toc
    
    % Decoding the variable sized blocks and placing them:
    tic;
    CCbrec = degroup_and_decode(Seq_acbr, dc_br, lbl_b, Fname, n/2, 2);
    toc
    
    % Resizing the Cb component:-
    % Subsampling ratio used: 4:2:0.
    tic;
    %CCbs = zeros(m/2,n/2);
    clear CCbr;
    k = 1;
    for i = 1:2:m
        l = 1;
        for j = 1:2:n
            CCbr(i:i+1,j:j+1) = CCbrec(k,l);
            l = l + 1;
        end
        k = k + 1;
    end
    %figure, imshow(uint8(CCbr));
    toc
    
    %%%%%%%%%%%%%%%%%%% Decoding the Cr Component %%%%%%%%%%%%%%%%%%%%%%
    
    %%% Decoding DC:
    % Decoding using Arithmetic QM coding;
    tic;
    Seq_dcrr = QMdecoder(Bout_dcr, length(Seq_dcr));
    toc
    
    % Reconstructing DC coefficients using Seq_dcr:
    tic;
    dc_rr = dc_seq_re1(Seq_dcrr, length(dc_r));
    toc
    
    %%% Decoding AC:
    % Decoding using Arithmetic QM coding;
    tic;
    Seq_acrr = QMdecoder(Bout_acr, length(Seq_acr));
    toc
    
   % Decoding the variable sized blocks and placing them:
    tic;
    CCrrec = degroup_and_decode(Seq_acrr, dc_rr, lbl_r, Fname, n/2, 2);
    toc
    
    % Resizing the Cb component:-
    % Subsampling ratio used: 4:2:0.
    tic;
    %CCbs = zeros(m/2,n/2);
    clear CCrr;
    k = 1;
    for i = 1:2:m
        l = 1;
        for j = 1:2:n
            CCrr(i:i+1,j:j+1) = CCrrec(k,l);
            l = l + 1;
        end
        k = k + 1;
    end
    %figure, imshow(uint8(CCrr));
    toc
    
    % YCbCr to RGB conversion:
    
    clear Cres;
    Cres(:,:,1) = CYrec;
    Cres(:,:,2) = CCbr;
    Cres(:,:,3) = CCrr;
    
    clear Crgb;
    %Cres = uint8(Cres);
    figure, imshow(uint8(Cres));
    %[Crgb, Cxyz] = Labtoxyzrgb(Cres);
    Crgb = lab2rgb(Cres,'OutputType','uint8');
    %Crgb = ycbcr2rgb(Cres);
    figure, imshow(Crgb);
    title('Reconstructed Image');
    
    OPBaseFileName = sprintf('OPImage%d.bmp', k10);
    OPFullFileName = fullfile(OPFolder,OPBaseFileName);
    imwrite(Crgb,OPFullFileName,'bmp');
    Tdecomp(k10) = toc(Tstrt1)
    
    %%
    tic;
    % Calculating MSE:-
    %load('Aoriginal.mat');
    clear Er;
    Er = (A - Crgb).^2;
    mse(k10) = (sum(sum(sum(Er))))/(m*n*3);
    psnr(k10) = 10*log10(65025/mse(k10));
    toc
    
    % Calculating SSIM:-
    tic;
    Ssim(k10) = 0;
    for k = 1:3
        Ssim(k10) = Ssim(k10) + ssim(A(:,:,k), Crgb(:,:,k));
    end
    Ssim(k10) = Ssim(k10)/3;
    
    tic;
    % Calculating CR:

    no_bits = (length(Bout_dcy)+length(Bout_acy)+length(Bout_dcb)+length(Bout_acb)+length(Bout_dcr)+length(Bout_acr))*8;
    
    %CR(k10) = (m*n*8*3)/no_bits;
    
    no_bits1 = length(Code_arr_y) + length(Code_arr_b) + length(Code_arr_r);
    
    CRe(k10) = (m*n*8*3)/(no_bits+no_bits1);
    
    % Saving the compressed bitstream arrays to a text file:-
    OPBFileName = fullfile(OPFolder,'bitstream.txt');
    fid = fopen(OPBFileName,'w+');
    fprintf(fid,'%s\r\n %s\r\n %s\r\n %s\r\n %s\r\n %s\r\n','Code_arr_y', Code_arr_y, 'Code_arr_b', Code_arr_b, 'Code_arr_r', Code_arr_r);
    fprintf(fid,'%s\r\n', 'Bout_dcy');
    fprintf(fid,'%6d', Bout_dcy);
    fprintf(fid,'\r\n');
    fprintf(fid,'%s\r\n', 'Bout_acy');
    fprintf(fid,'%6d', Bout_acy);
    fprintf(fid,'\r\n');
    fprintf(fid,'%s\r\n', 'Bout_dcb');
    fprintf(fid,'%6d', Bout_dcb);
    fprintf(fid,'\r\n');
    fprintf(fid,'%s\r\n', 'Bout_acb');
    fprintf(fid,'%6d', Bout_acb);
    fprintf(fid,'\r\n');
    fprintf(fid,'%s\r\n', 'Bout_dcr');
    fprintf(fid,'%6d', Bout_dcr);
    fprintf(fid,'\r\n');
    fprintf(fid,'%s\r\n', 'Bout_acr');
    fprintf(fid,'%6d', Bout_acr);
    fprintf(fid,'\r\n');  
    fclose(fid);
    
    
    % Compressed Size in kB:
    Siz(k10) = no_bits/(8*1024);
    toc
end
%xlswrite(OPXFileName,Szfull);
sheet=1;
xlswrite(OPXFileName,Siz',sheet); % Compressed size in kBs.
sheet=2;
xlswrite(OPXFileName,Tcomp',sheet); % Encoding time in seconds.
sheet=3;
xlswrite(OPXFileName,Tdecomp',sheet); % Decoding time in seconds.
sheet=4;
xlswrite(OPXFileName,CRe',sheet); % Compression Ratio
sheet=5;
xlswrite(OPXFileName,mse',sheet); % MSE
sheet=6;
xlswrite(OPXFileName,psnr',sheet); % PSNR
sheet=7;
xlswrite(OPXFileName,Ssim',sheet); %SSIM