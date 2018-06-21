function [gd] = ddmm2gd(ddmm)
%
% function [gd] = ddmm2gd(ddmm)
%
%   input ddmm : lati(ddmm), longi(ddmm), alti(m), 1 x 3 matrix
%   output gd = lati(degree), longi(degree), alti(meter)
%
%
%   coded by Joonseong Gim, Feb 26, 2016
la = fix(ddmm(1)/100) + (ddmm(1)/100-fix(ddmm(1)/100))*100/60;
lo = fix(ddmm(2)/100) + (ddmm(2)/100-fix(ddmm(2)/100))*100/60;
h = ddmm(3);
gd = [la lo h];