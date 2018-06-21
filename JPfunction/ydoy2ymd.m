function [day, month, year] = ydoy2ymd(yy, ddd)
%
%function [day, month, year] = ydoy2ymd(yy, ddd)
%
% DO: Convert YY & DOY to DAY, MONTH, YEAR
%
% <input>   yy: year in 2-digit
%           ddd: DOY(day of Year)
%
% <output>  day: Day in 2-digit
%           month: Month in 2-digit
%           year: Year in 4-digit
%
% eg. [d, m, y] = ydoy2ymd(03,100)
%
%
% Copyright: Kwan-Dong Park; February 8, 2009
%

if yy < 80
    year = 2000 + yy;
else
    year = 1900 + yy;
end

months = [1 2 3 4 5 6 7 8 9 10 11 12];

switch (year)
    case {1980, 1984, 1988, 1992, 1996, 2000, 2004, 2008, 2012, 2016, 2020}
        days = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
        ndays = 366;
    otherwise
        days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
        ndays = 365;
end

if ddd <= 0 || ddd > ndays
    return;
end

day = ddd;

i = 1;
while day > days(i)
    day = day - days(i);
    i = i + 1;
end
month = months(i);
