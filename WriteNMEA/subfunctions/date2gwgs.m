function [gw, gs] = date2gwgs(yy, mo, dd, h, m, s)
%
%function [gw, gs] = date2gwgs(yy, mo, dd, h, m, s)
%
% DO: Calculates [GW, GS] from Year, Month, Day, Hour, Minute, Seccond
%       - GW: GPS Week Number
%       - GS: GPS Week Seconds
%
% Copyright: Coded by Kwan-Dong Park, October 2001, Harvard-Smithsonian CfA
%
jd = Date2jd(yy, mo, dd, h, m, s);
[gw, gs] = jd2gwgs(jd);

