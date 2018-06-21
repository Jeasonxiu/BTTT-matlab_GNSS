function [jd, mjd] = ymd2jd(GPS_Time)
%
%Usage: ymd2jd(2005,9,28,11,23,10)
%inputs: year, month, day, hour, min,sec
%outputs: julian day, modified julian day
%
%by: Jihyun Ha
%      Kookimin Univ. 
%Creation Date: Oct 4, 2005
%

format long g;

year = GPS_Time(1);
month = GPS_Time(2);
day = GPS_Time(3);
hour = GPS_Time(4);
min = GPS_Time(5);
sec = GPS_Time(6);

if nargin==3
    hour=0;
    min=0;
    sec=0;
elseif nargin==4
    min=0;
    sec=0;
end

if month < 3
    year=year-1;
    month=month+12;
end

ut=hour+(min/60)+(sec/3600);
jd=fix(365.25*year)+fix(30.6001*(month+1))+day+(ut/24)+1720981.5;
mjd=jd-2400000.5;
