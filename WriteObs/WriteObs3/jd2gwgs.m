function [gw,gs] = jd2gwgs(JD)
%
% DO: Convert Julian Day (JD) to GPS Week and GPS Week Seconds
%
% <input>   JD: Julian Day
%
% <output>  gw: GPS Week
%           gs: GPS Week Seconds
%
% Copyright: Kwan-Dong Park, December 10, 2013 @LDEO
%            - original: [gps_time] by Kai Borre

a = floor(JD + .5);
b = a + 1537;
c = floor((b - 122.1) / 365.25);
e = floor(365.25*c);
f = floor((b - e) / 30.6001);
d = b - e - floor(30.6001*f) + rem(JD + .5, 1);
day_of_week = rem(floor(JD + .5), 7);
gw = floor((JD - 2444244.5) / 7);
gs = (rem(d, 1) + day_of_week + 1) * 86400;
gs = mod(gs, (86400 * 7));

