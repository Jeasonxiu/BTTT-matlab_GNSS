function [GSA] = GSAmat2(gpgsa)
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

GSA = zeros(length(gpgsa),12);
if ~isempty(gpgsa)
    for i = 1: length(gpgsa)
        line = cell2mat(gpgsa(i));
        nlength=length(line);
        [msg,pgn,prn1,prn2,prn3,prn4,prn5,prn6,prn7,prn8,prn9,prn10,prn11,prn12,dop1,dop2,dop3] =  ...
            strread(line(8:nlength-4),'%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','delimiter',',');
        GSA(i,:) = [prn1,prn2,prn3,prn4,prn5,prn6,prn7,prn8,prn9,prn10,prn11,prn12];
    end
else
    GSA =[];
end
prnlist = GSA(find(GSA(:,:) > 1));
GSA = sort(prnlist);

 
                