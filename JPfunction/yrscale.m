function [yrs] = yrscale(ddd, yy)
%function [yrscale] = yrscale(ddd, yy)
%
% Convert DDD & YY to floating Year
%
% Copyright: Kwan-Dong Park, February 29, 2008
%
leapYear = [1996 2000 2004 2008 2012 2016 2020];
n = length(leapYear);

for k = 1:n-1
    if yy >= leapYear(k) && yy < leapYear(k+1)
        lastLeapYear = leapYear(k);
        break
    end
end

fracDay = 1/365.25;

nYr = yy - lastLeapYear;
if nYr == 0
    nDay = ddd - 1;
else
    nDay = 366 + (nYr - 1) * 365 + ddd - 1;
end

yrs = lastLeapYear + nDay * fracDay;

str_yrs=num2str(yrs,'%4.4f');
num_yrs=str2num(str_yrs);
yrs=num_yrs; 