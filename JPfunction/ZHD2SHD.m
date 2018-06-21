function [SHD] = ZHD2SHD(gw, gs, pos, ele, ZHD)
%
%function [SHD] = ZHD2SHD(gw, gs, pos, ele)
%
% DO: Convert ZHD to SHD
%
% <input>   gw/gs: GPS Week Number and GPS Week Second
%           pos: Site position
%           ele: elevation angle in degrees
%           ZHD: ZHD in [mm]
%
% <output>  SHD: SHD in [m]
%
% Copyright: Kwan-Dong Park, 12/18/2014
%
[gmfh, dum] = GMF(gw, gs, pos, ele);
SHD = ZHD*gmfh;     %: [mm]
SHD = SHD/1000.;    %: Conversion to [m]