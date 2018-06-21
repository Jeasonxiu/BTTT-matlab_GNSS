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

%% Boehm의 GPT를 사용하기 위한 준비 작업
%* Convert Position to Lat, Lon, Hgt
gd = xyz2gd(pos); lat = gd(1); lon = gd(2); hgt = gd(3); % hgt for height: [m]
%* Convert GW/GS to MJD
mjd = gwgs2mjd(gw, gs);
%% GPT로 압력(P) 산출
[p, dum, dum] = gpt_v1(mjd, deg2rad(lat), deg2rad(lon), hgt); %: lat/lon[RAD], hgt[M]
%% 압력기반으로 ZHD 계산, 그리고 간단 사상함수 1/sin(el) 적용
% ZHD = (2.2779 * p) / (1 - 0.00266 * cos(lat) - 0.00028 * hgt);      % 원본 함수
ZHD = (2.2779 * p) / (1 - 0.00266 * cos(lat) - 0.00028 * hgt/1000);   % 수정 함수
delTrop = ZHD;
% delTrop = ZHD/sind(ele);      %: [km]
% delTrop = delTrop/1000.;      %: Conversion to [METER]    
% delTrop = ZHD/1000.;      %: Conversion to [METER]          % 원본 함수
