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
% 4/19/14: lat�� 2*lat�� ����, hgt�� 1000���� ������ ������ km�� ����
% 4/19/14: ���ܻ���Լ� 1/sin(el)�� �ƴ� GMF ����
% 12/18/14: ������� �������� �ʰ� ZHD�� ����

%% Boehm�� GPT�� ����ϱ� ���� �غ� �۾�
%* Convert Position to Lat, Lon, Hgt
gd = xyz2gd(pos); lat = gd(1); lon = gd(2); hgt = gd(3); % hgt for height: [m]
%* Convert GW/GS to MJD
mjd = gwgs2mjd(gw, gs);
%% GPT�� �з�(P) ����
[p, temp, undu] = gpt_v1(mjd, deg2rad(lat), deg2rad(lon), hgt); %: lat/lon[RAD], hgt[M]
%% �з±������ ZHD ���, �׸��� GMF ����Լ� 1/sin(el) ����
ZHD = (2.277 * p) / (1 - 0.00266 * cosd(2*lat) - 0.00028 * hgt/1000);
