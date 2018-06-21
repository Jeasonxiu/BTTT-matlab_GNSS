function [iers] = iers(jd);
    
%    input jd = Julian date
%    
%    output EOPC04 = [xp, yp, dUT1, LOD, dPsi, dEps]
%    reference : http://www.iers.org/IERS/EN/DataProducts/EarthOrientationData/eop.html

[yy,mo,dd,hh,mm,ss] = jd2date(jd);

load('iers.mat')
indexyear = find(iers(:,1) == yy);
iers = iers(indexyear,:);
indexmonth = find(iers(:,2) == mo);
iers = iers(indexmonth,:);
indexday = find(iers(:,3) == dd);
iers = iers(indexday,:);
iers = iers(1,5:10);