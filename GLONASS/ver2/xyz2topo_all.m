function topo = xyz2topo_all(xyzs,lat,lon)
%
% XYZ2TOPO: TRANSFORM FROM XYZ TO TOPOCENTRIC USING LATITUDE AND LONGITUDE
%
% xyz2topo(xyzs,lat,lon)
%   * lat:lon - given in degrees. they are converted into 'radians' in this subroutine.
%

% Conversion to Radians
lat=lat*pi/180.;
lon=lon*pi/180.;

% Cosine and Sines of Latitude and Longitude
cos_lat=cos(lat);
sin_lat=sin(lat);
cos_lon=cos(lon);
sin_lon=sin(lon);

% Conversion to Topocentric 
%   - Vallado's Fundamentals of Astrodynamics and Applications
% 수정전
% xx=xyzs(1); yy=xyzs(2); zz=xyzs(3);
% 
% topo(1,1)=-sin_lat*cos_lon*xx - sin_lat*sin_lon*yy + cos_lat*zz; % n
% topo(1,2)=        -sin_lon*xx +         cos_lon*yy;              % e
% topo(1,3)= cos_lat*cos_lon*xx + cos_lat*sin_lon*yy + sin_lat*zz; % v

%  수정
%  모든 xyzs에 대해서 nev변환이 가능하도록 수정
xx=xyzs(:,1); yy=xyzs(:,2); zz=xyzs(:,3);

for i = 1:length(xx)

topo(i,1)=-sin_lat*cos_lon*xx(i) - sin_lat*sin_lon*yy(i) + cos_lat*zz(i); % n
topo(i,2)=        -sin_lon*xx(i) +         cos_lon*yy(i);                 % e
topo(i,3)= cos_lat*cos_lon*xx(i) + cos_lat*sin_lon*yy(i) + sin_lat*zz(i); % v

end