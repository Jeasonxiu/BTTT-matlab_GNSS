function[x,y,z] = azel2xyz(la,lo,ht)
%
%function[x,y,z] = azel2xyz(la,lo,ht)
%
% <input>   la, lo, ht : 위경도, 고도
%
% <output>  x, y, z : xyz 좌표
%
% Copyright: Jinyi Kim, November 26, 2014@INHA University
%
format long g;

D2R = pi/180;
la = la*D2R;
lo = lo*D2R;

a = 6378137.0;
b = 6356752.314;
N = a^2/sqrt((a^2)*(cos(la)^2)+(b^2)*(sin(la)^2));

x = (N+ht)*cos(la)*cos(lo);
y = (N+ht)*cos(la)*sin(lo);
z = ((b^2/a^2)*N+ht)*sin(la);
