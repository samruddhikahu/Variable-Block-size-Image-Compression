function [blk,idx] = Arithmetic_Decoding_blk_csf_lab_lin_ksii(Seq, FileName, dc, blksz, flag)
% load('Nmat.mat');
% load('Cnmat.mat');
%load('Table.mat');
% if(flag==1)
%     mat = Nmat;
% else
%     mat = Cnmat;
% end
m = blksz(1,1);
n = blksz(1,2);

% CSF based quantization matrices for different block sizes:-
if(flag==1)
    if(m==8)&&(n==8)
        load(fullfile(FileName,'quant_lin8.mat'));
        %quant8(1,1) = 20;
        quant = quant_lin;
    elseif(m==8)&&(n==16)
        load(fullfile(FileName,'quant_lin816.mat'));
        %quant8_16(1,1) = 25;
        quant = quant_lin;
    elseif(m==8)&&(n==24)
        load(fullfile(FileName,'quant_lin824.mat'));
        %quant8_16(1,1) = 25;
        quant = quant_lin;
    elseif(m==8)&&(n==32)
        load(fullfile(FileName,'quant_lin832.mat'));
        %quant8_16(1,1) = 25;
        quant = quant_lin;
    elseif(m==16)&&(n==8)
        load(fullfile(FileName,'quant_lin168.mat'));
        %quant16_8(1,1) = 25;
        quant = quant_lin;
    elseif(m==16)&&(n==16)
        load(fullfile(FileName,'quant_lin16.mat'));
        %quant16(1,1) = 30;
        quant = quant_lin;
    elseif(m==16)&&(n==24)
        load(fullfile(FileName,'quant_lin1624.mat'));
        %quant16_24(1,1) =35;
        quant = quant_lin;
    elseif(m==16)&&(n==32)
        load(fullfile(FileName,'quant_lin1632.mat'));
        %quant16_24(1,1) =35;
        quant = quant_lin;
    elseif(m==24)&&(n==8)
        load(fullfile(FileName,'quant_lin248.mat'));
        %quant24_16(1,1) = 35;
        quant = quant_lin;
    elseif(m==24)&&(n==16)
        load(fullfile(FileName,'quant_lin2416.mat'));
        %quant24_16(1,1) = 35;
        quant = quant_lin;
    elseif(m==24)&&(n==24)
        load(fullfile(FileName,'quant_lin24.mat'));
        %quant24(1,1) = 40;
        quant = quant_lin;
    elseif(m==24)&&(n==32)
        load(fullfile(FileName,'quant_lin2432.mat'));
        %quant24_32(1,1) = 45;
        quant = quant_lin;
    elseif(m==32)&&(n==8)
        load(fullfile(FileName,'quant_lin328.mat'));
        %quant32_24(1,1) = 45;
        quant = quant_lin;
    elseif(m==32)&&(n==16)
        load(fullfile(FileName,'quant_lin3216.mat'));
        %quant32_24(1,1) = 45;
        quant = quant_lin;
    elseif(m==32)&&(n==24)
        load(fullfile(FileName,'quant_lin3224.mat'));
        %quant32_24(1,1) = 45;
        quant = quant_lin;
    elseif(m==32)&&(n==32)
        load(fullfile(FileName,'quant_lin32.mat'));
        %quant32(1,1) = 50;
        quant = quant_lin;
    end
    %quant = quant*2;
else
    if(m==8)&&(n==8)
        load(fullfile(FileName,'quant_lin_ab8.mat'));
        %quant8(1,1) = 20;
        quant = quant_lin;
    elseif(m==8)&&(n==16)
        load(fullfile(FileName,'quant_lin_ab816.mat'));
        %quant8_16(1,1) = 25;
        quant = quant_lin;
    elseif(m==8)&&(n==24)
        load(fullfile(FileName,'quant_lin_ab824.mat'));
        %quant8_16(1,1) = 25;
        quant = quant_lin;
    elseif(m==8)&&(n==32)
        load(fullfile(FileName,'quant_lin_ab832.mat'));
        %quant8_16(1,1) = 25;
        quant = quant_lin;
    elseif(m==16)&&(n==8)
        load(fullfile(FileName,'quant_lin_ab168.mat'));
        %quant16_8(1,1) = 25;
        quant = quant_lin;
    elseif(m==16)&&(n==16)
        load(fullfile(FileName,'quant_lin_ab16.mat'));
        %quant16(1,1) = 30;
        quant = quant_lin;
    elseif(m==16)&&(n==24)
        load(fullfile(FileName,'quant_lin_ab1624.mat'));
        %quant16_24(1,1) =35;
        quant = quant_lin;
    elseif(m==16)&&(n==32)
        load(fullfile(FileName,'quant_lin_ab1632.mat'));
        %quant16_24(1,1) =35;
        quant = quant_lin;
    elseif(m==24)&&(n==8)
        load(fullfile(FileName,'quant_lin_ab248.mat'));
        %quant24_16(1,1) = 35;
        quant = quant_lin;
    elseif(m==24)&&(n==16)
        load(fullfile(FileName,'quant_lin_ab2416.mat'));
        %quant24_16(1,1) = 35;
        quant = quant_lin;
    elseif(m==24)&&(n==24)
        load(fullfile(FileName,'quant_lin_ab24.mat'));
        %quant24(1,1) = 40;
        quant = quant_lin;
    elseif(m==24)&&(n==32)
        load(fullfile(FileName,'quant_lin_ab2432.mat'));
        %quant24_32(1,1) = 45;
        quant = quant_lin;
    elseif(m==32)&&(n==8)
        load(fullfile(FileName,'quant_lin_ab328.mat'));
        %quant32_24(1,1) = 45;
        quant = quant_lin;
    elseif(m==32)&&(n==16)
        load(fullfile(FileName,'quant_lin_ab3216.mat'));
        %quant32_24(1,1) = 45;
        quant = quant_lin;
    elseif(m==32)&&(n==24)
        load(fullfile(FileName,'quant_lin_ab3224.mat'));
        %quant32_24(1,1) = 45;
        quant = quant_lin;
    elseif(m==32)&&(n==32)
        load(fullfile(FileName,'quant_lin_ab32.mat'));
        %quant32(1,1) = 50;
        quant = quant_lin;
    end
    %quant = quant*1.5;
end
%quant = 2*quant;
% Nmat_new = round(Nmat_new);

% Arithmetic Decoding of AC Coefficients:
p = 0;
k = 1;
while(k<=(m*n - 1))
    p = p + 1;
    D = str2num(Seq(p));
    if(D==1)
        Ir(1,k:(m*n - 1)) = 0;
        break;
    else
        p = p + 1;
        D = str2num(Seq(p));
        while(D==0)
            Ir(1,k) = 0;
            k = k + 1;
            p = p + 1;
            D = str2num(Seq(p));
        end
        % Decode V:
        % Decode Sign of V:
        p = p + 1;
        D = str2num(Seq(p));
        if(D==1)
            Sgn = 'neg';
        else
            Sgn = 'pos';
        end
        % Decode log2(Sz):
        n1 = 0;
        p = p + 1;
        D = str2num(Seq(p));
        while(D)
            n1 = n1 + 1;
            p = p + 1;
            D = str2num(Seq(p));
        end
        n1 = n1 - 1;
        if(n1==0)
            Sz = 1;
        elseif(n1>0)
            %base = 2^n;
            Szr = bin2dec(Seq(p+1:p+n1));
            Sz = Szr + (2^n1);
            p = p + n1;
        else
            Sz = 0;
        end
        abV = Sz + 1;
        if(Sgn=='pos')
            V = abV;
        else
            V = -abV;
        end
        
        Ir(1,k) = V;
        k = k + 1;
    end
end
idx = p;

Irev_y = [dc Ir];

% Zigzag reordering:
sm = m+n-1;
s = 1;
bq(1,1) = Irev_y(s,1);
z = 2;
k = 1;
if(m>n)
    n1 = min(m,n);
    m1 = max(m,n);
    for y = 3:sm
        if(y<(n1+2))
            if(mod(y,2)==1)
                k = 1;
                while((y-k)>0)
                    bq(k,y-k) = Irev_y(s,z);
                    k = k + 1;
                    z = z + 1;
                end
            else
                l = 1;
                while(((y-l)>0)&&(y<=n1))
                    bq(y-l,l) = Irev_y(s,z);
                    l = l + 1;
                    z = z + 1;
                end
            end
        elseif((y>=(n1+2))&&(y<(m1+2)))
            if(mod(y,2)==1)
                l = n;
                k = y-l;
                while((y-k)>0)
                    bq(k,y-k) = Irev_y(s,z);
                    k = k + 1;
                    z = z + 1;
                end
            else
                l = 1;
                while(((y-l)>=0)&&(l<=n1))
                    bq(y-l,l) = Irev_y(s,z);
                    l = l + 1;
                    z = z + 1;
                end
            end
        else
            if(mod(y,2)==1)
                l = n;
                k = y-l;
                while(((y-k)>0)&&(k<=m1))
                    bq(k,y-k) = Irev_y(s,z);
                    k = k + 1;
                    z = z + 1;
                end
            else
                l = m;
                k = y-l;
                while(((y-k)>=0)&&(k<=n1))
                    bq(y-k,k) = Irev_y(s,z);
                    k = k + 1;
                    z = z + 1;
                end
            end
        end
    end
elseif(m<n)
    n1 = min(m,n);
    m1 = max(m,n);
    for y = 3:sm
        if(y<(n1+2))
            if(mod(y,2)==1)
                k = 1;
                while((y-k)>0)
                    bq(k,y-k) = Irev_y(s,z);
                    k = k + 1;
                    z = z + 1;
                end
            else
                l = 1;
                while(((y-l)>0)&&(y<=n1))
                    bq(y-l,l) = Irev_y(s,z);
                    l = l + 1;
                    z = z + 1;
                end
            end
        elseif((y>=(n1+2))&&(y<(m1+2)))
            if(mod(y,2)==1)
                k = 1;
                while(((y-k)>0)&&(k<=n1))
                    bq(k,y-k) = Irev_y(s,z);
                    k = k + 1;
                    z = z + 1;
                end
            else
                k = n1;
                l = y-k;
                while((y-l)>0)
                    bq(y-l,l) = Irev_y(s,z);
                    l = l + 1;
                    z = z + 1;
                end
            end
        else
            if(mod(y,2)==1)
                l = n;
                k = y-l;
                while(((y-k)>0)&&(k<=n1))
                    bq(k,y-k) = Irev_y(s,z);
                    k = k + 1;
                    z = z + 1;
                end
            else
                l = m;
                k = y-l;
                while(((y-k)>=0)&&(k<=m1))
                    bq(y-k,k) = Irev_y(s,z);
                    k = k + 1;
                    z = z + 1;
                end
            end
        end
    end
elseif(m==n)
    for y = 3:sm
        if(y<(m+2))
            if(mod(y,2)==1)
                k = 1;
                while((y-k)>0)
                    bq(k,y-k) = Irev_y(s,z);
                    k = k + 1;
                    z = z + 1;
                end
            else
                l = 1;
                while((y-l)>0)
                    bq(y-l,l) = Irev_y(s,z);
                    l = l + 1;
                    z = z + 1;
                end
            end
        else
            if(mod(y,2)==1)
                l = m;
                k = y-l;
                while((l-k)>=0)
                    bq(k,y-k) = Irev_y(s,z);
                    k = k + 1;
                    z = z + 1;
                end
            else
                k = n;
                l = y-k;
                while((k-l)>=0)
                    bq(y-l,l) = Irev_y(s,z);
                    l = l + 1;
                    z = z + 1;
                end
            end
        end
    end
end          
bq(m,n) = Irev_y(s,z);

if(flag==1)
    sh = 128;
else
    sh = 0;
end

% Reverse Quantization:
bdct = bq; %.*quant;
% Inverse DCT:
bshifted = idct2(bdct);
% Level shifting:
blk = bshifted + sh;
end
