function xyz = gd2xyz(gd)
%
%function [xyz] = gd2xyz2(gd)
%
% DO: Convert GD(geodetic) coordinates into XYZ
%
% <input>   gd: Geodetic cooridnates (lat/lon/ht) 1X3 row vector
%
% <output>  xyz: XYZ coordinates (X/Y/Z) 1X3 row vector
%
% Coded by: Kwan-Dong Park, April 07, 2007
% --- Modifications ---
% 3/29/14: 출력을 row vector로 변경
%

%% 변수 [la : latitude; lo : longitude; h  : height]
toRad = pi/180;
la = gd(1); la = la * toRad;
lo = gd(2); lo = lo * toRad;
h  = gd(3);
%% WGS 84, GRS 80 제원: a와 f
% GRS80: a = 6378137.0, 1/f = 298.257222101
% WGS84: a = 6378137.0, 1/f = 298.257223563 
a = 6378137.0; f = 1/298.257223563; 
%% 좌표 변환 과정
e = sqrt( 2*f - f^2);
N = a / sqrt(1 - e^2 * sin(la)^2);
b = a * (1 - f);

x = (N + h) * cos(la) * cos(lo);
y = (N + h) * cos(la) * sin(lo);
z = (b^2/a^2 * N + h) * sin(la);

xyz = [x y z];
