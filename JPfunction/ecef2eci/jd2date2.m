function [yy,mo,dd,hh,mm,ss] = jd2date2(JD)
% 
%function [yy,mo,dd,hh,mm,ss] = jd2date(JD)
%
% DO: Convert Julian Day to DATE(Calendar Day)
%
% <input>   JD: Julian Day
%
% <output>  yy, mo, dd: year/month/day
%           hh, mm, ss: hour/minute/second
%
% Copyright: Kwan-Dong Park, December 10, 2013 @LDEO
%           - original by Jihyun Ha, 2005
% Modyfied by Joonseong,
%

a = floor(JD + .5);
b = a + 1537;
c = floor((b - 122.1) / 365.25);
d = floor(365.25*c);
e = floor((b - d) / 30.6001);

dd = b - d - floor(30.6001*e) + rem(JD + 0.5, 1);
hh = rem(dd, 1)*24;
mm = rem(hh, 1)*60;
ss = rem(mm, 1)*60;

% ss = round(ss);
ss = ss;
mm = floor(mm);
hh = floor(hh);
dd = floor(dd);
mo = e - 1 - 12*floor(e/14);
yy = c - 4715 - floor((7 + mo)/10);
