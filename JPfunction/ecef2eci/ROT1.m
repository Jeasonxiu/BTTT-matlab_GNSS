function [ROT1] = ROT1(X)
    
%    input X = Degree
    
ROT1 = [1, 0, 0; 0, cosd(X), sind(X); 0, -sind(X), cosd(X)];
