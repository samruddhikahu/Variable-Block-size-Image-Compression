function Fname = quantization_matrix_adaptive_ksii(const)
% clc;
% clear all;
% close all;

Fname = sprintf('Lab Quantization Matrices IIQ Linear%0.2d',const);
QFolder = fullfile(Fname);
if ~exist(QFolder, 'dir')
    mkdir(QFolder)
end

mat_no = [8, 816, 824, 832, 168, 16, 1624, 1632, 248, 2416, 24, 2432, 328, 3216, 3224, 32];
mat_size = [8,8;
            8,16;
            8,24;
            8,32;
            16,8;
            16,16;
            16,24;
            16,32;
            24,8;
            24,16;
            24,24;
            24,32;
            32,8;
            32,16;
            32,24;
            32,32];
% For luminance quantization matrices:-
base = 'coeff_max_L';
%DC_val_L = [10, 10, 10, 15, 15, 15, 20, 20, 20, 25];
DC_val = [8, 11, 14, 16, 11, 16, 20, 23, 14, 20, 24, 28, 16, 23, 28, 32];

for i = 1:16
    clear cname quant_lin;
    cname = [base num2str(mat_no(i))];
    load(cname);
    quant_lin = quantization_matrix_csf_adap(mat_size(i,1),mat_size(i,2),50,coeff_max,const);
    quant_lin(1,1) = DC_val(1,i);
    quant_lin = round(quant_lin);
    BFName = sprintf('quant_lin%d.mat',mat_no(i));
    QBFileName = fullfile(QFolder,BFName);
    save(QBFileName,'quant_lin');
end

% For chrominance quantization matrices:-
base = 'coeff_max_ab';
%DC_val_ab = [5, 7, 7, 10, 13, 13, 15, 18, 18, 20];

for i = 1:16
    clear cname quant_lin;
    cname = [base num2str(mat_no(i))];
    load(cname);
    quant_lin = quantization_matrix_csf_adap(mat_size(i,1),mat_size(i,2),100,coeff_max,const);
    quant_lin(1,1) = DC_val(1,i);
    quant_lin = round(quant_lin);
    BFName = sprintf('quant_lin_ab%d.mat',mat_no(i));
    QBFileName = fullfile(QFolder,BFName);
    save(QBFileName,'quant_lin');
end
end