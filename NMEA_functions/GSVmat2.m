function [GSV_sort] = GSVmat2(gpgsv)
%
% function [GSV] = GSVmat2(gpgsv)
%
%   Read the gpgsv per epoch, get GSV matrix from gpgsv
%
%   input gpgsv : gpgsv/epoch
%
%   Example : [] = GSVmat2(gpgsv)
%
%   coded by Joonseong Gim, DEC 12, 2016

GSV = zeros(1,16);
if ~isempty(gpgsv)
    k=1;
    for i = 1: length(gpgsv)
        line = cell2mat(gpgsv(i));
        nlength=length(line);
        if line(nlength-3) == ','
            add= strcat(',', line(nlength-3:nlength));
            line(nlength-3:nlength+1) = add;
            nlength=length(line);
        end
        index = findstr(line,',');
        if length(index) >= 19
            [pgn,sq,nos,prn1,el1,az1,snr1,prn2,el2,az2,snr2,prn3,el3,az3,snr3,prn4,el4,az4,snr4] =  ...
                strread(line(8:nlength-3),'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','delimiter',',');
            GSV(1,k:k+15) = [prn1,el1,az1,snr1,prn2,el2,az2,snr2,prn3,el3,az3,snr3,prn4,el4,az4,snr4];
            k= k+16;
        elseif length(index) < 19 & length(index) >= 14
            [pgn,sq,nos,prn1,el1,az1,snr1,prn2,el2,az2,snr2,prn3,el3,az3,snr3] =  ...
                strread(line(8:nlength-3),'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','delimiter',',');
            GSV(1,k:k+15) = [prn1,el1,az1,snr1,prn2,el2,az2,snr2,prn3,el3,az3,snr3,0,0,0,0];
            k= k+16;
        elseif length(index) < 14 & length(index) >= 11
            [pgn,sq,nos,prn1,el1,az1,snr1,prn2,el2,az2,snr2] =  ...
                strread(line(8:nlength-3),'%f%f%f%f%f%f%f%f%f%f%f','delimiter',',');
            GSV(1,k:k+15) = [prn1,el1,az1,snr1,prn2,el2,az2,snr2,0,0,0,0,0,0,0,0];
            k= k+16;
        elseif length(index) < 11
            [pgn,sq,nos,prn1,el1,az1,snr1] =  ...
                strread(line(8:nlength-3),'%f%f%f%f%f%f%f','delimiter',',');
            GSV(1,k:k+15) = [prn1,el1,az1,snr1,0,0,0,0,0,0,0,0,0,0,0,0];
            k= k+16;
            
        end
    end
    %     GSV=GSV_sort;
else
    GSV_sort = [];
end
if ~isempty(GSV)
    i=1; k=1;
    for i =1:nos
        GSV_sort(i,:) = GSV(k:k+3);
        k = k + 4;
    end
end






