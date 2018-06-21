function [hh,mm,ss,x,y,z,la,lo,qi,nSats,ht] = NEWreadGGA3(line2) 
%
%function [hh,mm,ss,la,lo,ht] = readGGA(line2)
%
% <input>   line2: �ش� NMEA �������� ��
%
% <output>  hh : ��
%           mm : ��
%           ss : ��
%           x  : ������ǥ�� X ��ǥ
%           y  : ������ǥ�� Y ��ǥ
%           z  : ������ǥ�� Z ��ǥ
%
% Copyright: __________, November 26, 2014@INHA University
% ���� : ���浵 ���� �� �Ҽ��� 4�ڸ��� �д� ���� �ذ�, ������(Oct 18 , 2016)
% line2= GGA;
format long g;
line2 = line2(1:length(line2)-3);
nlength = length(line2);

[utc,lat,ns,lon,ew,qi,nSats,hdop,htGeoid,unit1,dgeoid,unit2,update,ref_id] =  ...
   strread(line2(8:nlength),'%s%f%s%f%s%d%d%f%f%s%f%s%f%d','delimiter',','); 


% [utc,lat,ns,lon,ew,qi,nSats,hdop,htGeoid,unit1,dgeoid,unit2] =  ...
%    strread(line2(8:nlength),'%f%f%s%f%s%d%d%f%f%s%f%s','delimiter',','); 

% GGA 1 : UTC
% if utc < 1e5
%     utcStr = num2str(utc);
%     hh = str2num(utcStr(1));
%     mm = str2num(utcStr(2:3));
%     ss = str2num(utcStr(4:end));
% else
%     utcStr = num2str(utc);
%     hh = str2num(utcStr(1:2));
%     mm = str2num(utcStr(3:4));
%     ss = str2num(utcStr(5:end));
% end
utcStr = cell2mat(utc);
hh = str2num(utcStr(1:2));
mm = str2num(utcStr(3:4));
ss = str2num(utcStr(5:end));
% GGA 2 : Latitute
latStr = num2str(lat,'%15.8f');
degLat = str2num(latStr(1:2));
minLat = str2num(latStr(3:end));
la = degLat + minLat/60;
% GGA 4 : Longitude
lonStr = num2str(lon,'%15.8f');
degLon = str2num(lonStr(1:3));
minLon = str2num(lonStr(4:end));
lo = degLon + minLon/60;
% GGA 9 & 11 : Height
if ~isempty(dgeoid)
    ht = htGeoid+dgeoid;
else
    ht = htGeoid;
end
% ������ǥ�� -> ������ǥ�� ��ȯ

[x,y,z] = azel2xyz(la,lo,ht);
% fprintf(' %s %s %15.9f %15.9f %8.9f\n',latStr, lonStr, la, lo, ht);
