function quant = quantization_matrix_csf_adap(M,N,range,coeff_max,const)
% This function is used to calculate quantization matrix for different
% color spaces
quant = zeros(M,N);
Nn = sqrt(M*N);
fmax = ceil(30*sqrt((M-1)^2 + (N-1)^2)/(Nn*1.5));
for m = 1:M
    for n = 1:N
        f = 30*sqrt((m-1)^2 + (n-1)^2)/(Nn*1.5);  % Assuming pixel size is 1.5 min
        if(f==0)
            CSF = 1;
        else
            CSF = const*(f-fmax);%100*sqrt(f)*exp(-0.13*f);
        end
        if(m==1)||(n==1)
            Norm = 2^(-2.5); 
            %Norm = 1;
            thresh = 1/(Norm*CSF);
        else
            if(m<=n)
                OTF = exp(-9.5*((m-1)/(n-1))^2);
            elseif(n<=m)
                OTF = exp(-9.5*((n-1)/(m-1))^2);
            end
            if(OTF<0.5)
                OTF = 0.5;
            end
            Norm = 2^(-2);
            %Norm = 1; %OTF = 1;
            thresh = 1/(OTF*Norm*CSF);
        end
        quant(m,n) = min(thresh*range,coeff_max(m,n));  
    end
end
end