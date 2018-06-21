function [gw,gs]=ymd2gpstime(year, month, day, hour, min, sec)
%
%Usage: ymd2gpstime(2005,9,28,11,23,10)
%inputs: year, month, day, hour, min,sec
%outputs: GPS Week Number[num], GPS Day Number[num], GPS Week Second[num], Day of Year[str]
%
%by: Jihyun Ha
%      Kookimin Univ. 
%Creation Date: Oct 4, 2005
%

format long g;

[jd, mjd] = ymd2jd([year,month,day,hour,min,sec]);

gw=fix((jd-2444244.5)/7);
if hour <12 
    gd=round((((jd-2444244.5)/7)-gw)*7);
elseif hour>=12
    gd=fix((((jd-2444244.5)/7)-gw)*7);
end
gs=round((gd*86400)+(hour*3600)+(min*60)+sec);

switch (year)
case {1980, 1984, 1988,1992,1996,2000,2004,2008, 2012}
    days=[31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
    ndays=366;
otherwise
    days=[31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
    ndays=365;
end

months=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];
doy=day;

i=1;
while month ~= months(i)
    doy=doy+days(i);
    i=i+1;
    if i>12
        break;
    end
end

if doy<10
    doy=strcat('00',num2str(doy));
elseif doy<100
    doy=strcat('0', num2str(doy));
end