function [gw,gd] = doy2gwgd(doy,year)

% input data ¼³Á¤
% year = 2014;
% doy = 25;

switch (year)
    case {1980,1984,1988,1992,1996,2000,2004,2008,2012,2016,2020,2024}
        days = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
    otherwise
        days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
end

tdays = 0;
for i = 1: 12
    tdays = tdays + days(i);
    if doy - tdays < 0
        break;
    end
end
month = i;
day = doy-sum(days(1:i-1));

[JD] = date2jd(year, month, day, 0, 0, 0);
gd = rem(floor(JD + .5), 7);
if gd == 6
    gd = 0;
else gd = gd + 1;
end

[gw, gs] = date2gwgs(year, month, day, 0, 0, 0);