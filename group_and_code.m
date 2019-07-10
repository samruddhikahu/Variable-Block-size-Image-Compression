function [Seq_acy, dc_y, lbl_y, lbl_arr_y] = group_and_code(CY, Idx, Fname, flag)
[m,n] = size(CY);
blksz = zeros(m/8,n/8);
lbl_y = zeros(m/8,n/8);
k = 1;
i = 1;
Seq_acy = [];
lbl_arr_y = [];
p = 1;

while(all(all(lbl_y))==0)
    l = 1;
    j = 1;
    while(j<=(n/8))
        clear Seq;
        r = min(find(lbl_y(:,j)==0));
        if(isempty(r))
            %r = 65;
            %flag = 'rc';
            j = j + 1;
            continue;
            %break;
        end
        i = r;
        if(i<=((m/8)-3))&&(j<=((n/8)-3))&&(any(any(lbl_y(i:i+3,j:j+3)))==0)&&(Idx(i,j)==Idx(i+1,j))&&(Idx(i,j)==Idx(i,j+1))&&(Idx(i,j)==Idx(i+1,j+1))&&(Idx(i,j)==Idx(i+2,j))&&(Idx(i,j)==Idx(i,j+2))&&(Idx(i,j)==Idx(i+1,j+2))...
                &&(Idx(i,j)==Idx(i+2,j+1))&&(Idx(i,j)==Idx(i+2,j+2))&&(Idx(i,j)==Idx(i+3,j))&&(Idx(i,j)==Idx(i+3,j+1))&&(Idx(i,j)==Idx(i+3,j+2))&&(Idx(i,j)==Idx(i,j+3))...
                &&(Idx(i,j)==Idx(i+1,j+3))&&(Idx(i,j)==Idx(i+2,j+3))&&(Idx(i,j)==Idx(i+3,j+3))
            blksz(k,l,1:2) = [32 32];
            lbl_y(i:i+3,j:j+3) = 16;
            lbl_arr_y(p) = 16;
            a1 = 8*(i-1) + 1;
            b1 = 8*(i-1) + 32;
            a2 = 8*(j-1) + 1;
            b2 = 8*(j-1) + 32;
            blk = CY(a1:b1,a2:b2);
            [Seq, dc] = Arithmetic_Coding_blk_csf_lab_lin_ksii(blk,Fname,flag);
            
            j = j + 4;
        elseif(i<=((m/8)-3))&&(j<=((n/8)-2))&&(any(any(lbl_y(i:i+3,j:j+2)))==0)&&(Idx(i,j)==Idx(i+1,j))&&(Idx(i,j)==Idx(i,j+1))&&(Idx(i,j)==Idx(i+1,j+1))&&(Idx(i,j)==Idx(i+2,j))&&(Idx(i,j)==Idx(i,j+2))&&(Idx(i,j)==Idx(i+1,j+2))...
                &&(Idx(i,j)==Idx(i+2,j+1))&&(Idx(i,j)==Idx(i+2,j+2))&&(Idx(i,j)==Idx(i+3,j))&&(Idx(i,j)==Idx(i+3,j+1))&&(Idx(i,j)==Idx(i+3,j+2))
            blksz(k,l,1:2) = [32 24];
            lbl_y(i:i+3,j:j+2) = 15;
            lbl_arr_y(p) = 15;
            a1 = 8*(i-1) + 1;
            b1 = 8*(i-1) + 32;
            a2 = 8*(j-1) + 1;
            b2 = 8*(j-1) + 24;
            blk = CY(a1:b1,a2:b2);
            [Seq, dc] = Arithmetic_Coding_blk_csf_lab_lin_ksii(blk,Fname,flag);
            
            j = j + 3;
        elseif(i<=((m/8)-2))&&(j<=((n/8)-3))&&(any(any(lbl_y(i:i+2,j:j+3)))==0)&&(Idx(i,j)==Idx(i+1,j))&&(Idx(i,j)==Idx(i,j+1))&&(Idx(i,j)==Idx(i+1,j+1))&&(Idx(i,j)==Idx(i+2,j))&&(Idx(i,j)==Idx(i,j+2))&&(Idx(i,j)==Idx(i+1,j+2))...
                &&(Idx(i,j)==Idx(i+2,j+1))&&(Idx(i,j)==Idx(i+2,j+2))&&(Idx(i,j)==Idx(i,j+3))&&(Idx(i,j)==Idx(i+1,j+3))&&(Idx(i,j)==Idx(i+2,j+3))
            blksz(k,l,1:2) = [24 32];
            lbl_y(i:i+2,j:j+3) = 14;
            lbl_arr_y(p) = 14;
            a1 = 8*(i-1) + 1;
            b1 = 8*(i-1) + 24;
            a2 = 8*(j-1) + 1;
            b2 = 8*(j-1) + 32;
            blk = CY(a1:b1,a2:b2);
            [Seq, dc] = Arithmetic_Coding_blk_csf_lab_lin_ksii(blk,Fname,flag);
            
            j = j + 4;
        elseif(i<=((m/8)-2))&&(j<=((n/8)-2))&&(any(any(lbl_y(i:i+2,j:j+2)))==0)&&(Idx(i,j)==Idx(i+1,j))&&(Idx(i,j)==Idx(i,j+1))&&(Idx(i,j)==Idx(i+1,j+1))&&(Idx(i,j)==Idx(i+2,j))&&(Idx(i,j)==Idx(i,j+2))&&(Idx(i,j)==Idx(i+1,j+2))...
                &&(Idx(i,j)==Idx(i+2,j+1))&&(Idx(i,j)==Idx(i+2,j+2))
            blksz(k,l,1:2) = [24 24];
            lbl_y(i:i+2,j:j+2) = 13;
            lbl_arr_y(p) = 13;
            a1 = 8*(i-1) + 1;
            b1 = 8*(i-1) + 24;
            a2 = 8*(j-1) + 1;
            b2 = 8*(j-1) + 24;
            blk = CY(a1:b1,a2:b2);
            [Seq, dc] = Arithmetic_Coding_blk_csf_lab_lin_ksii(blk,Fname,flag);
            
            j = j + 3;
        elseif(i<=((m/8)-3))&&(j<=((n/8)-1))&&(any(any(lbl_y(i:i+3,j:j+1)))==0)&&(Idx(i,j)==Idx(i+1,j))&&(Idx(i,j)==Idx(i,j+1))&&(Idx(i,j)==Idx(i+1,j+1))&&(Idx(i,j)==Idx(i+2,j))...
                &&(Idx(i,j)==Idx(i+2,j+1))&&(Idx(i,j)==Idx(i+3,j))&&(Idx(i,j)==Idx(i+3,j+1))
            blksz(k,l,1:2) = [32 16];
            lbl_y(i:i+3,j:j+1) = 12;
            lbl_arr_y(p) = 12;
            a1 = 8*(i-1) + 1;
            b1 = 8*(i-1) + 32;
            a2 = 8*(j-1) + 1;
            b2 = 8*(j-1) + 16;
            blk = CY(a1:b1,a2:b2);
            [Seq, dc] = Arithmetic_Coding_blk_csf_lab_lin_ksii(blk,Fname,flag);
            
            j = j + 2;
        elseif(i<=((m/8)-1))&&(j<=((n/8)-3))&&(any(any(lbl_y(i:i+1,j:j+3)))==0)&&(Idx(i,j)==Idx(i+1,j))&&(Idx(i,j)==Idx(i,j+1))&&(Idx(i,j)==Idx(i+1,j+1))...
                &&(Idx(i,j)==Idx(i,j+2))&&(Idx(i,j)==Idx(i,j+3))&&(Idx(i,j)==Idx(i+1,j+2))&&(Idx(i,j)==Idx(i+1,j+3))
            blksz(k,l,1:2) = [16 32];
            lbl_y(i:i+1,j:j+3) = 11;
            lbl_arr_y(p) = 11;
            a1 = 8*(i-1) + 1;
            b1 = 8*(i-1) + 16;
            a2 = 8*(j-1) + 1;
            b2 = 8*(j-1) + 32;
            blk = CY(a1:b1,a2:b2);
            [Seq, dc] = Arithmetic_Coding_blk_csf_lab_lin_ksii(blk,Fname,flag);
            
            j = j + 4;
        elseif(i<=((m/8)-2))&&(j<=((n/8)-1))&&(any(any(lbl_y(i:i+2,j:j+1)))==0)&&(Idx(i,j)==Idx(i+1,j))&&(Idx(i,j)==Idx(i,j+1))&&(Idx(i,j)==Idx(i+1,j+1))&&(Idx(i,j)==Idx(i+2,j))&&(Idx(i,j)==Idx(i+2,j+1))
            blksz(k,l,1:2) = [24 16];
            lbl_y(i:i+2,j:j+1) = 10;
            lbl_arr_y(p) = 10;
            a1 = 8*(i-1) + 1;
            b1 = 8*(i-1) + 24;
            a2 = 8*(j-1) + 1;
            b2 = 8*(j-1) + 16;
            blk = CY(a1:b1,a2:b2);
            [Seq, dc] = Arithmetic_Coding_blk_csf_lab_lin_ksii(blk,Fname,flag);
            
            j = j + 2;
        elseif(i<=((m/8)-1))&&(j<=((n/8)-2))&&(any(any(lbl_y(i:i+1,j:j+2)))==0)&&(Idx(i,j)==Idx(i+1,j))&&(Idx(i,j)==Idx(i,j+1))&&(Idx(i,j)==Idx(i+1,j+1))&&(Idx(i,j)==Idx(i,j+2))&&(Idx(i,j)==Idx(i+1,j+2))
            blksz(k,l,1:2) = [16 24];
            lbl_y(i:i+1,j:j+2) = 9;
            lbl_arr_y(p) = 9;
            a1 = 8*(i-1) + 1;
            b1 = 8*(i-1) + 16;
            a2 = 8*(j-1) + 1;
            b2 = 8*(j-1) + 24;
            blk = CY(a1:b1,a2:b2);
            [Seq, dc] = Arithmetic_Coding_blk_csf_lab_lin_ksii(blk,Fname,flag);
            
            j = j + 3;
        elseif(i<=((m/8)-3))&&(any(any(lbl_y(i:i+3,j)))==0)&&(Idx(i,j)==Idx(i+1,j))&&(Idx(i,j)==Idx(i+2,j))&&(Idx(i,j)==Idx(i+3,j))
            blksz(k,l,1:2) = [32 8];
            lbl_y(i:i+3,j) = 8;
            lbl_arr_y(p) = 8;
            a1 = 8*(i-1) + 1;
            b1 = 8*(i-1) + 32;
            a2 = 8*(j-1) + 1;
            b2 = 8*(j-1) + 8;
            blk = CY(a1:b1,a2:b2);
            [Seq, dc] = Arithmetic_Coding_blk_csf_lab_lin_ksii(blk,Fname,flag);
            
            j = j + 1;
        elseif(i<=((m/8)-1))&&(j<=((n/8)-1))&&(any(any(lbl_y(i:i+1,j:j+1)))==0)&&(Idx(i,j)==Idx(i+1,j))&&(Idx(i,j)==Idx(i,j+1))&&(Idx(i,j)==Idx(i+1,j+1))
            blksz(k,l,1:2) = [16 16];
            lbl_y(i:i+1,j:j+1) = 7;
            lbl_arr_y(p) = 7;
            a1 = 8*(i-1) + 1;
            b1 = 8*(i-1) + 16;
            a2 = 8*(j-1) + 1;
            b2 = 8*(j-1) + 16;
            blk = CY(a1:b1,a2:b2);
            [Seq, dc] = Arithmetic_Coding_blk_csf_lab_lin_ksii(blk,Fname,flag);
            
            j = j + 2;
        elseif(j<=((n/8)-3))&&(any(any(lbl_y(i,j:j+3)))==0)&&(Idx(i,j)==Idx(i,j+1))&&(Idx(i,j)==Idx(i,j+2))&&(Idx(i,j)==Idx(i,j+3))
            blksz(k,l,1:2) = [8 32];
            lbl_y(i,j:j+3) = 6;
            lbl_arr_y(p) = 6;
            a1 = 8*(i-1) + 1;
            b1 = 8*(i-1) + 8;
            a2 = 8*(j-1) + 1;
            b2 = 8*(j-1) + 32;
            blk = CY(a1:b1,a2:b2);
            [Seq, dc] = Arithmetic_Coding_blk_csf_lab_lin_ksii(blk,Fname,flag);
            
            j = j + 4;
        elseif(i<=((m/8)-2))&&(any(any(lbl_y(i:i+2,j)))==0)&&(Idx(i,j)==Idx(i+1,j))&&(Idx(i,j)==Idx(i+2,j))
            blksz(k,l,1:2) = [24 8];
            lbl_y(i:i+2,j) = 5;
            lbl_arr_y(p) = 5;
            a1 = 8*(i-1) + 1;
            b1 = 8*(i-1) + 24;
            a2 = 8*(j-1) + 1;
            b2 = 8*(j-1) + 8;
            blk = CY(a1:b1,a2:b2);
            [Seq, dc] = Arithmetic_Coding_blk_csf_lab_lin_ksii(blk,Fname,flag);
            
            j = j + 1;
        elseif(j<=((n/8)-2))&&(any(any(lbl_y(i,j:j+2)))==0)&&(Idx(i,j)==Idx(i,j+1))&&(Idx(i,j)==Idx(i,j+2))
            blksz(k,l,1:2) = [8 24];
            lbl_y(i,j:j+2) = 4;
            lbl_arr_y(p) = 4;
            a1 = 8*(i-1) + 1;
            b1 = 8*(i-1) + 8;
            a2 = 8*(j-1) + 1;
            b2 = 8*(j-1) + 24;
            blk = CY(a1:b1,a2:b2);
            [Seq, dc] = Arithmetic_Coding_blk_csf_lab_lin_ksii(blk,Fname,flag);
            
            j = j + 3;
        elseif(i<=((m/8)-1))&&(any(any(lbl_y(i:i+1,j)))==0)&&(Idx(i,j)==Idx(i+1,j))
            blksz(k,l,1:2) = [16 8];
            lbl_y(i:i+1,j) = 3;
            lbl_arr_y(p) = 3;
            a1 = 8*(i-1) + 1;
            b1 = 8*(i-1) + 16;
            a2 = 8*(j-1) + 1;
            b2 = 8*(j-1) + 8;
            blk = CY(a1:b1,a2:b2);
            [Seq, dc] = Arithmetic_Coding_blk_csf_lab_lin_ksii(blk,Fname,flag);
            
            j = j + 1;
        elseif(j<=((n/8)-1))&&(any(any(lbl_y(i,j:j+1)))==0)&&(Idx(i,j)==Idx(i,j+1))
            blksz(k,l,1:2) = [8 16];
            lbl_y(i,j:j+1) = 2;
            lbl_arr_y(p) = 2;
            a1 = 8*(i-1) + 1;
            b1 = 8*(i-1) + 8;
            a2 = 8*(j-1) + 1;
            b2 = 8*(j-1) + 16;
            blk = CY(a1:b1,a2:b2);
            [Seq, dc] = Arithmetic_Coding_blk_csf_lab_lin_ksii(blk,Fname,flag);
            
            j = j + 2;
        else
            blksz(k,l,1:2) = [8 8];
            lbl_y(i,j) = 1;
            lbl_arr_y(p) = 1;
            a1 = 8*(i-1) + 1;
            b1 = 8*(i-1) + 8;
            a2 = 8*(j-1) + 1;
            b2 = 8*(j-1) + 8;
            blk = CY(a1:b1,a2:b2);
            [Seq, dc] = Arithmetic_Coding_blk_csf_lab_lin_ksii(blk,Fname,flag);
            
            j = j + 1;
        end
        Seq_acy = [Seq_acy Seq];
        %Seq_acy{p} = Seq;
        
        % Grouping the DC coefficients:-
        dc_y(p) = dc;
        
        p = p + 1;
        l = l + 1;
    end
    %j = 1;
    k = k + 1;
    %fl = 1;
    i = i + 1;
end

end