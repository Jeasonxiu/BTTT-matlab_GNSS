function [az,el] = xyz2azel(xyz,lat,lon)
%
%function [az,el] = xyz2azel(xyz,lat,lon)
%
%   Convert (dX, dY, dZ) to (Azimuth, Elevation)
%
%       - Coded by Kwan-Dong Park; November 2011; Inha University
%       - Input: XYZ = 3 X 1 column vector
%                lat/lon = Degress   
%

%* convert column vector to row vector
[a, b] = size(xyz);
if a==3
    xyz = xyz';
end
topo = xyz2topo(xyz,lat,lon);

N = topo(1);
E = topo(2);
V = topo(3);

proj = sqrt(N^2 + E^2);
el = atan2(V, proj) * 180/pi;
az = angle(N + i * E) * 180/pi;

if az < 0
    az = 360 + az;
end