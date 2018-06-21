function [delTrop] = TropGPTh(pos, gw, gs, ele)
%function [delTrop] = TropGPTh(pos, gw, gs, ele)
%
% DO: Get Hydrostatic Delay based on GPT
%
% <input>   pos: site position in XYZ [m]
%           gw/gs: GPS Week number and seconds, used for MJD computation
%           ele: elevation angle in degrees
%
% <output>  delTrop: Zenith hydrostatic Delay [m]
%
% Copyright: Kwan-Dong Park, December 2013 @LDEO
%

%% Boehm�� GPT�� ����ϱ� ���� �غ� �۾�
%* Convert Position to Lat, Lon, Hgt
gd = xyz2gd(pos); lat = gd(1); lon = gd(2); hgt = gd(3); % hgt for height: [m]
%* Convert GW/GS to MJD
mjd = gwgs2mjd(gw, gs);
%% GPT�� �з�(P) ����
[p, dum, dum] = gpt_v1(mjd, deg2rad(lat), deg2rad(lon), hgt); %: lat/lon[RAD], hgt[M]
%% �з±������ ZHD ���, �׸��� ���� ����Լ� 1/sin(el) ����
% ZHD = (2.2779 * p) / (1 - 0.00266 * cos(lat) - 0.00028 * hgt);      % ���� �Լ�
ZHD = (2.2779 * p) / (1 - 0.00266 * cos(lat) - 0.00028 * hgt/1000);   % ���� �Լ�
delTrop = ZHD;
% delTrop = ZHD/sind(ele);      %: [km]
% delTrop = delTrop/1000.;      %: Conversion to [METER]    
% delTrop = ZHD/1000.;      %: Conversion to [METER]          % ���� �Լ�
