function [topo] = xyz2topo(xyzs, lat, lon)
%
%function [topo] = xyz2topo(xyzs, lat, lon)
%
% DO: Transform XYZ to topocentric NEV using latitude and longitude
%
%       xyzs: 3-column matrix
%       lat, lon: [degrees] - to be converted to radians in the routine
%
% Copyright: Kwan-Dong Park, revised on January 28, 2008
%            Modified on Dec 06, 2009
%               - to handle array input
%------------------------------------------------------------------------------

%% '도(degrees)'를 라디안으로 변경
lat = lat * pi / 180;
lon = lon * pi / 180;
%% 위경도의 cosine, sine
cos_lat = cos(lat);
sin_lat = sin(lat);
cos_lon = cos(lon);
sin_lon = sin(lon);
%% 좌표 변환, 참고 Vallado's Fundamentals of Astrodynamics and Applications
xx = xyzs(:,1); yy = xyzs(:,2); zz = xyzs(:,3);
topo(:,1) = -sin_lat*cos_lon.*xx - sin_lat*sin_lon.*yy + cos_lat.*zz;
topo(:,2) =         -sin_lon.*xx +         cos_lon.*yy;
% topo(:,3) =  cos_lat*cos_lon.*xx + cos_lat*sin_lon.*yy + sin_lat.*zz;
topo(:,3) =  cos_lat*cos_lon.*xx + cos_lat*sin_lon.*yy + sin_lat.*zz;
