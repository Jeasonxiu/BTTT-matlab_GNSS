function [gmfh, gmfw] = GMF(gw, gs, pos, el)
%
%function [m_h, m_w] = GMF(gw, gs, pos, el)
%
% DO: Derive GMF hydrostatic and wet mapping function
%       - GMF: Global Mapping Function
%
% <input>   gw/gs: GPS Week Number and Weekk Seconds
%           pos: site position in ECEF [m]
%           el: elevation angle [DEG]
%
% <output>  gmfh/gmfw: GMF Hydrostatic/Wet mapping function
%
% Copyright: Kwan-Dong Park, January 23, 2014 @LDEO
%
% 기본적으로 Boehm이 만든 gmf_f_hu를 구동하기 위한 입출력 기능에 지나지 않음
%

%% MJD 환산
dmjd = gwgs2mjd(gw, gs);
%% 위경도-고도 환산(단위 rad/rad/m)
llh = xyz2gd(pos);
dlat = deg2rad(llh(1));
dlon = deg2rad(llh(2));
dhgt = llh(3);
%% Zenith Distance 환산(단위 rad)
zd = deg2rad(90 - el);
%% gmf_f_hu 호출해서 사상함수 결정
[gmfh, gmfw] = gmf_f_hu (dmjd, dlat, dlon, dhgt, zd);
