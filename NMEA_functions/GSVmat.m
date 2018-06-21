function [GSV] = GSVmat(gpgsv)
%
% function [GSV] = GSVmat(gpgsv)
%
%   Read the gpgsv per epoch, get GSV matrix from gpgsv
%   
%   input gpgsv : gpgsv/epoch
%
%   Example : [] = GSVmat(gpgsv)
%
%   coded by Joonseong Gim, Jan 25, 2016

GSV = zeros(length(gpgsv),19);
for i = 1: length(gpgsv)
    line = cell2mat(gpgsv(i));
    index = findstr(line,',');
    for j = 2:length(index)
        if index(j) - index(j-1) == 1
            GSV(i,j-1) = 0;
        else
            GSV(i,j-1) = str2num(line(index(j-1)+1:index(j)-1));
            if j==length(index) 
                if line(length(line)-3:length(line)) ~= ','
                    GSV(i,j) = str2num(line(index(j)+1:index(j)+2));
                else
                    GSV(i,j) = 0;
                end
            end
        end
    end
    
end


 
                