function [hh,mm,ss,x,y,z,la,lo,ht,qi] = NEWreadGGA(line2) 
%
%function [hh,mm,ss,la,lo,ht] = readGGA(line2)
%
% <input>   line2: ÇØ´ç NMEA µ¥ÀÌÅÍÀÇ ÁÙ
%
% <output>  hh : ½Ã
%           mm : ºÐ
%           ss : ÃÊ
%           x  : Á÷°¢ÁÂÇ¥°è X ÁÂÇ¥
%           y  : Á÷°¢ÁÂÇ¥°è Y ÁÂÇ¥
%           z  : Á÷°¢ÁÂÇ¥°è Z ÁÂÇ¥
%
% Copyright: __________, November 26, 2014@INHA University
%

nlength = length(line2);

[utc,lat,ns,lon,ew,qi,nSats,hdop,htGeoid,unit1,dgeoid,unit2,update,ref_id] =  ...
   strread(line2(8:nlength),'%f%f%s%f%s%d%d%f%f%s%f%s%f%d','delimiter',','); 

% GGA 1 : UTC
if utc < 1e5
    utcStr = num2str(utc);
    hh = str2num(utcStr(1));
    mm = str2num(utcStr(2:3));
    ss = str2num(utcStr(4:end));
else
    utcStr = num2str(utc);
    hh = str2num(utcStr(1:2));
    mm = str2num(utcStr(3:4));
    ss = str2num(utcStr(5:end));
end
% GGA 2 : Latitute
lat = lat*1e8;
latStr = num2str(lat);
degLat = str2num(latStr(1:2));
minLat = str2num(latStr(3:end));
minLat = minLat/1e8;
la = degLat + minLat/60;
% GGA 4 : Longitude
lon = lon*1e8;
lonStr = num2str(lon);
degLon = str2num(lonStr(1:3));
minLon = str2num(lonStr(4:end));
minLon = minLon/1e8;
lo = degLon + minLon/60;
% GGA 9 & 11 : Height
ht = htGeoid+dgeoid;
% ÁöÇüÁÂÇ¥°è -> Á÷°¢ÁÂÇ¥°è º¯È¯
[x,y,z] = azel2xyz(la,lo,ht);
