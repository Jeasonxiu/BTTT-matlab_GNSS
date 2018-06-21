function [d, m, s] = deg2dms(deg)

d = floor(deg);
m = floor((deg-d)*60);
s = rem((deg-d)*60,1) * 60;

fprintf('%3d %2d %11.8f\n',d,m,s)