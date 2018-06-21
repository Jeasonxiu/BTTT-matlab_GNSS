function [ZHD] = TropGPTh(pos, gw, gs)
%function [ZHD] = TropGPTh(pos, gw, gs, ele)
%
% DO: Get Zenith Hydrostatic Delay based on GPT
%
% <input>   pos: site position in XYZ [m]
%           gw/gs: GPS Week number and seconds, used for MJD computation
%      
% <output>  ZHD: Zenith hydrostatic Delay [mm]
%
% Copyright: Kwan-Dong Park, December 2013 @LDEO
% --- Modifications ---
% 4/19/14: lat를 2*lat로 변경, hgt를 1000으로 나눠서 단위를 km로 변경
% 4/19/14: 간단사상함수 1/sin(el)이 아닌 GMF 적용
% 12/18/14: 매핑펑션 적용하지 않고 ZHD만 추출

%% Boehm의 GPT를 사용하기 위한 준비 작업
%* Convert Position to Lat, Lon, Hgt
gd = xyz2gd(pos); lat = gd(1); lon = gd(2); hgt = gd(3); % hgt for height: [m]
%* Convert GW/GS to MJD
mjd = gwgs2mjd(gw, gs);
%% GPT로 압력(P) 산출
[p, temp, undu] = gpt_v1(mjd, deg2rad(lat), deg2rad(lon), hgt); %: lat/lon[RAD], hgt[M]
%% 압력기반으로 ZHD 계산, 그리고 GMF 사상함수 1/sin(el) 적용
ZHD = (2.277 * p) / (1 - 0.00266 * cosd(2*lat) - 0.00028 * hgt/1000);
