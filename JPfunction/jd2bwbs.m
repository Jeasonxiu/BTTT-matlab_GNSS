function [bw,bs] = jd2bwbs(JD)
%
% DO: Convert Julian Day (JD) to BDS Week and BDS Week Seconds
%
% <input>   JD: Julian Day
%
% <output>  bw: BDS Week
%           bs: BDS Seconds of Week (SOW)
%
% Copyright: Kwan-Dong Park, 1/24/2015
%          

a = floor(JD + .5);
b = a + 1537;
c = floor((b - 122.1) / 365.25);
e = floor(365.25*c);
f = floor((b - e) / 30.6001);
d = b - e - floor(30.6001*f) + rem(JD + .5, 1);
day_of_week = rem(floor(JD + .5), 7);
bw = floor((JD - 2453736.5) / 7); %: 위성항법시스템에 따라 변경!
bs = (rem(d, 1) + day_of_week + 1) * 86400;
bs = mod(bs, (86400 * 7));

