function [dotR] = dotR(X)
    
%    input X = Degree
    
dotR = [-sind(X), -cosd(X), 0; cosd(X), -sind(X), 0;, 0, 0, 0];