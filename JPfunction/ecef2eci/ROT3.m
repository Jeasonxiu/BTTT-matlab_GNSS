function [ROT3] = ROT3(X)
    
%    input X = Degree
    
ROT3 = [cosd(X), sind(X), 0; -sind(X), cosd(X), 0; 0, 0, 1];
