function CYrec = degroup_and_decode(Seq_acyr, dc_yr, lbl_y, Fname, n, flag)
%lbl_y1 = lbl_y;
    i = 1; p = 1;
    idx = 1;
    %clear blksz;
    while(any(any(lbl_y))==1)
        l = 1; j = 1;
        while(j<=(n/8))
            r = min(find(lbl_y(:,j)));
            if(isempty(r))
                j = j + 1;
                continue;
            end
            i = r;
            clear blk k;            
            if(lbl_y(i,j)==1)
                blksz = [8 8];
                dc = dc_yr(p);
                [blk,k] = Arithmetic_Decoding_blk_csf_lab_lin_ksii(Seq_acyr(idx:end), Fname, dc, blksz, flag);
                a1 = 8*(i-1) + 1;
                b1 = 8*(i-1) + 8;
                a2 = 8*(j-1) + 1;
                b2 = 8*(j-1) + 8;
                CYrec(a1:b1,a2:b2) = blk;
                lbl_y(i,j) = 0;
                j = j + 1;
            elseif(lbl_y(i,j)==2)
                blksz = [8 16];
                dc = dc_yr(p);
                [blk,k] = Arithmetic_Decoding_blk_csf_lab_lin_ksii(Seq_acyr(idx:end), Fname, dc, blksz, flag);
                a1 = 8*(i-1) + 1;
                b1 = 8*(i-1) + 8;
                a2 = 8*(j-1) + 1;
                b2 = 8*(j-1) + 16;
                CYrec(a1:b1,a2:b2) = blk;
                lbl_y(i,j:j+1) = 0;
                j = j + 2;
            elseif(lbl_y(i,j)==3)
                blksz = [16 8];
                dc = dc_yr(p);
                [blk,k] = Arithmetic_Decoding_blk_csf_lab_lin_ksii(Seq_acyr(idx:end), Fname, dc, blksz, flag);
                a1 = 8*(i-1) + 1;
                b1 = 8*(i-1) + 16;
                a2 = 8*(j-1) + 1;
                b2 = 8*(j-1) + 8;
                CYrec(a1:b1,a2:b2) = blk;
                lbl_y(i:i+1,j) = 0;
                j = j + 1;
            elseif(lbl_y(i,j)==4)
                blksz = [8 24];
                dc = dc_yr(p);
                [blk,k] = Arithmetic_Decoding_blk_csf_lab_lin_ksii(Seq_acyr(idx:end), Fname, dc, blksz, flag);
                a1 = 8*(i-1) + 1;
                b1 = 8*(i-1) + 8;
                a2 = 8*(j-1) + 1;
                b2 = 8*(j-1) + 24;
                CYrec(a1:b1,a2:b2) = blk;
                lbl_y(i,j:j+2) = 0;
                j = j + 3;
            elseif(lbl_y(i,j)==5)
                blksz = [24 8];
                dc = dc_yr(p);
                [blk,k] = Arithmetic_Decoding_blk_csf_lab_lin_ksii(Seq_acyr(idx:end), Fname, dc, blksz, flag);
                a1 = 8*(i-1) + 1;
                b1 = 8*(i-1) + 24;
                a2 = 8*(j-1) + 1;
                b2 = 8*(j-1) + 8;
                CYrec(a1:b1,a2:b2) = blk;
                lbl_y(i:i+2,j) = 0;
                j = j + 1;
            elseif(lbl_y(i,j)==6)
                blksz = [8 32];
                dc = dc_yr(p);
                [blk,k] = Arithmetic_Decoding_blk_csf_lab_lin_ksii(Seq_acyr(idx:end), Fname, dc, blksz, flag);
                a1 = 8*(i-1) + 1;
                b1 = 8*(i-1) + 8;
                a2 = 8*(j-1) + 1;
                b2 = 8*(j-1) + 32;
                CYrec(a1:b1,a2:b2) = blk;
                lbl_y(i,j:j+3) = 0;
                j = j + 4;
            elseif(lbl_y(i,j)==7)
                blksz = [16 16];
                dc = dc_yr(p);
                [blk,k] = Arithmetic_Decoding_blk_csf_lab_lin_ksii(Seq_acyr(idx:end), Fname, dc, blksz, flag);
                a1 = 8*(i-1) + 1;
                b1 = 8*(i-1) + 16;
                a2 = 8*(j-1) + 1;
                b2 = 8*(j-1) + 16;
                CYrec(a1:b1,a2:b2) = blk;
                lbl_y(i:i+1,j:j+1) = 0;
                j = j + 2;
            elseif(lbl_y(i,j)==8)
                blksz = [32 8];
                dc = dc_yr(p);
                [blk,k] = Arithmetic_Decoding_blk_csf_lab_lin_ksii(Seq_acyr(idx:end), Fname, dc, blksz, flag);
                a1 = 8*(i-1) + 1;
                b1 = 8*(i-1) + 32;
                a2 = 8*(j-1) + 1;
                b2 = 8*(j-1) + 8;
                CYrec(a1:b1,a2:b2) = blk;
                lbl_y(i:i+3,j) = 0;
                j = j + 1;
            elseif(lbl_y(i,j)==9)
                blksz = [16 24];
                dc = dc_yr(p);
                [blk,k] = Arithmetic_Decoding_blk_csf_lab_lin_ksii(Seq_acyr(idx:end), Fname, dc, blksz, flag);
                a1 = 8*(i-1) + 1;
                b1 = 8*(i-1) + 16;
                a2 = 8*(j-1) + 1;
                b2 = 8*(j-1) + 24;
                CYrec(a1:b1,a2:b2) = blk;
                lbl_y(i:i+1,j:j+2) = 0;
                j = j + 3;
            elseif(lbl_y(i,j)==10)
                blksz = [24 16];
                dc = dc_yr(p);
                [blk,k] = Arithmetic_Decoding_blk_csf_lab_lin_ksii(Seq_acyr(idx:end), Fname, dc, blksz, flag);
                a1 = 8*(i-1) + 1;
                b1 = 8*(i-1) + 24;
                a2 = 8*(j-1) + 1;
                b2 = 8*(j-1) + 16;
                CYrec(a1:b1,a2:b2) = blk;
                lbl_y(i:i+2,j:j+1) = 0;
                j = j + 2;
            elseif(lbl_y(i,j)==11)
                blksz = [16 32];
                dc = dc_yr(p);
                [blk,k] = Arithmetic_Decoding_blk_csf_lab_lin_ksii(Seq_acyr(idx:end), Fname, dc, blksz, flag);
                a1 = 8*(i-1) + 1;
                b1 = 8*(i-1) + 16;
                a2 = 8*(j-1) + 1;
                b2 = 8*(j-1) + 32;
                CYrec(a1:b1,a2:b2) = blk;
                lbl_y(i:i+1,j:j+3) = 0;
                j = j + 4;
            elseif(lbl_y(i,j)==12)
                blksz = [32 16];
                dc = dc_yr(p);
                [blk,k] = Arithmetic_Decoding_blk_csf_lab_lin_ksii(Seq_acyr(idx:end), Fname, dc, blksz, flag);
                a1 = 8*(i-1) + 1;
                b1 = 8*(i-1) + 32;
                a2 = 8*(j-1) + 1;
                b2 = 8*(j-1) + 16;
                CYrec(a1:b1,a2:b2) = blk;
                lbl_y(i:i+3,j:j+1) = 0;
                j = j + 2;
            elseif(lbl_y(i,j)==13)
                blksz = [24 24];
                dc = dc_yr(p);
                [blk,k] = Arithmetic_Decoding_blk_csf_lab_lin_ksii(Seq_acyr(idx:end), Fname, dc, blksz, flag);
                a1 = 8*(i-1) + 1;
                b1 = 8*(i-1) + 24;
                a2 = 8*(j-1) + 1;
                b2 = 8*(j-1) + 24;
                CYrec(a1:b1,a2:b2) = blk;
                lbl_y(i:i+2,j:j+2) = 0;
                j = j + 3;
            elseif(lbl_y(i,j)==14)
                blksz = [24 32];
                dc = dc_yr(p);
                [blk,k] = Arithmetic_Decoding_blk_csf_lab_lin_ksii(Seq_acyr(idx:end), Fname, dc, blksz, flag);
                a1 = 8*(i-1) + 1;
                b1 = 8*(i-1) + 24;
                a2 = 8*(j-1) + 1;
                b2 = 8*(j-1) + 32;
                CYrec(a1:b1,a2:b2) = blk;
                lbl_y(i:i+2,j:j+3) = 0;
                j = j + 4;
            elseif(lbl_y(i,j)==15)
                blksz = [32 24];
                dc = dc_yr(p);
                [blk,k] = Arithmetic_Decoding_blk_csf_lab_lin_ksii(Seq_acyr(idx:end), Fname, dc, blksz, flag);
                a1 = 8*(i-1) + 1;
                b1 = 8*(i-1) + 32;
                a2 = 8*(j-1) + 1;
                b2 = 8*(j-1) + 24;
                CYrec(a1:b1,a2:b2) = blk;
                lbl_y(i:i+3,j:j+2) = 0;
                j = j + 3;
            elseif(lbl_y(i,j)==16)
                blksz = [32 32];
                dc = dc_yr(p);
                [blk,k] = Arithmetic_Decoding_blk_csf_lab_lin_ksii(Seq_acyr(idx:end), Fname, dc, blksz, flag);
                a1 = 8*(i-1) + 1;
                b1 = 8*(i-1) + 32;
                a2 = 8*(j-1) + 1;
                b2 = 8*(j-1) + 32;
                CYrec(a1:b1,a2:b2) = blk;
                lbl_y(i:i+3,j:j+3) = 0;
                j = j + 4;
            end
            p = p + 1;
            idx = idx + k;
            %hold on
            %imshow(uint8(CYrec));
        end
    end
    %figure, imshow(uint8(CYrec));

end