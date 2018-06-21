function xyz = topo2xyz(n,e,u,lat,lon)
%
%  xyz = topo2xyz(n,e,u,lat,lon)
%  
%  NEU ������ �������� XYZ ������ ��ȯ
%
%  input 
%       n, e, u : meter
%       lat, lon : ����, �浵 (degree ����)
%  output
%       xyz = [x;y;z] (���� meter)
%
% 2014.08.08 ������

% Cosine and Sines of Latitude and Longitude
cos_lat=cosd(lat);
sin_lat=sind(lat);
cos_lon=cosd(lon);
sin_lon=sind(lon);

T=[-sin_lat*cos_lon   -sin_lat*sin_lon   cos_lat;
   -sin_lon            cos_lon               0  ;
    cos_lat*cos_lon    cos_lat*sin_lon   sin_lat];

xyz=inv(T)*[n;e;u];


