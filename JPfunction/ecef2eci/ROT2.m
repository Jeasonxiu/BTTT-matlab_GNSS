function [ROT2] = ROT2(X)
    
%    input X = Degree
    
ROT2 = [cosd(X), 0, -sind(X); 0, 1, 0; sind(X), 0, cosd(X)];
