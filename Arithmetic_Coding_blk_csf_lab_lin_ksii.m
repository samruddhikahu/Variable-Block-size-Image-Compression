function [Seq, DC] = Arithmetic_Coding_blk_csf_lab_lin_ksii(Blk,FileName,flag)
% load('Nmat.mat');
% load('Cnmat.mat');
% if(flag==1)
%     mat = Nmat;
% else
%     mat = Cnmat;
% end
[m, n] = size(Blk);
sz = [m n];
quant = [];
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
%Nmat_new = round(Nmat_new);
if(flag==1)
    sh = 128;
else
    sh = 0;
end

b = double(Blk);
% Level shifting by -128
bshifted = b-sh;
% Calculating 2D-DCT:
bdct = dct2(bshifted);
% Quantizing using the normalization matrix defined above:
%bq = round(bdct./quant);
bq = round(bdct);
% Zigzag ordering:
sm = m+n-1;
s = 1;
I(s,1) = bq(1,1);
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
                    I(s,z) = bq(k,y-k);
                    k = k + 1;
                    z = z + 1;
                end
            else
                l = 1;
                while(((y-l)>0)&&(y<=n1))
                    I(s,z) = bq(y-l,l);
                    l = l + 1;
                    z = z + 1;
                end
            end
        elseif((y>=(n1+2))&&(y<(m1+2)))
            if(mod(y,2)==1)
                l = n;
                k = y-l;
                while((y-k)>0)
                    I(s,z) = bq(k,y-k);
                    k = k + 1;
                    z = z + 1;
                end
            else
                l = 1;
                while(((y-l)>=0)&&(l<=n1))
                    I(s,z) = bq(y-l,l);
                    l = l + 1;
                    z = z + 1;
                end
            end
        else
            if(mod(y,2)==1)
                l = n;
                k = y-l;
                while(((y-k)>0)&&(k<=m1))
                    I(s,z) = bq(k,y-k);
                    k = k + 1;
                    z = z + 1;
                end
            else
                l = m;
                k = y-l;
                while(((y-k)>=0)&&(k<=n1))
                    I(s,z) = bq(y-k,k);
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
                    I(s,z) = bq(k,y-k);
                    k = k + 1;
                    z = z + 1;
                end
            else
                l = 1;
                while(((y-l)>0)&&(y<=n1))
                    I(s,z) = bq(y-l,l);
                    l = l + 1;
                    z = z + 1;
                end
            end
        elseif((y>=(n1+2))&&(y<(m1+2)))
            if(mod(y,2)==1)
                k = 1;
                while(((y-k)>0)&&(k<=n1))
                    I(s,z) = bq(k,y-k);
                    k = k + 1;
                    z = z + 1;
                end
            else
                k = n1;
                l = y-k;
                while((y-l)>0)
                    I(s,z) = bq(y-l,l);
                    l = l + 1;
                    z = z + 1;
                end
            end
        else
            if(mod(y,2)==1)
                l = n;
                k = y-l;
                while(((y-k)>0)&&(k<=n1))
                    I(s,z) = bq(k,y-k);
                    k = k + 1;
                    z = z + 1;
                end
            else
                l = m;
                k = y-l;
                while(((y-k)>=0)&&(k<=m1))
                    I(s,z) = bq(y-k,k);
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
                    I(s,z) = bq(k,y-k);
                    k = k + 1;
                    z = z + 1;
                end
            else
                l = 1;
                while((y-l)>0)
                    I(s,z) = bq(y-l,l);
                    l = l + 1;
                    z = z + 1;
                end
            end
        else
            if(mod(y,2)==1)
                l = m;
                k = y-l;
                while((l-k)>=0)
                    I(s,z) = bq(k,y-k);
                    k = k + 1;
                    z = z + 1;
                end
            else
                k = n;
                l = y-k;
                while((k-l)>=0)
                    I(s,z) = bq(y-l,l);
                    l = l + 1;
                    z = z + 1;
                end
            end
        end
    end
end            
I(s,z) = bq(m,n);

% DC Coefficient:
DC = I(1,1);
% Arithmetic QM Coding:
Seq = [];
count = 0;
z = 1;
for i = 2:(m*n)
    if(I(1,i)==0)
        count = count + 1;
        last(z) = i;
        z = z + 1;
    else
        count = 0;
        z = 1;
    end
end

If = I(1,:);
if(count>0)
    %If(last(1)) = 10;
    if(last(1)<=(m*n))
        If(last(1):(m*n)) = [];
    end
end
Ig = If(:,2:end);
len = length(Ig);

k = 1;
while(k<=(m*n - 1))
    if(k>len)
        % Code 1:
        Seq = [Seq '1'];
        break;
    else
        % Code 0:
        Seq = [Seq '0'];
        
        V = Ig(k);
        while(V==0)
            % Code 0:
            Seq = [Seq '0'];
            k = k + 1;
            V = Ig(k);
        end
        
        % Code 1:
        Seq = [Seq '1'];
        % Encode V:
        if(V<0)
            Seq = [Seq '1'];
        else
            Seq = [Seq '0'];
        end
        Sz = abs(V) - 1;
        n1 = floor(log2(Sz));
        for i = 1:n1+1
            Seq = [Seq '1'];
        end
        Seq = [Seq '0'];
        if(n1>0)
            Szr = Sz - (2^n1);
            Seq = [Seq dec2bin(Szr,n1)];
        end
    end
    k = k + 1;
end

end
%toc