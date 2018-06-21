function [jd] = gwgs2jd(gw, gs)
% 
%   coded by Joonseong Gim, April 6, 2016
%   
[yy, mo, dd, hh, mm, ss] = gwgs2date(gw, gs);
[jd] = date2jd(yy, mo, dd, hh, mm, ss);
