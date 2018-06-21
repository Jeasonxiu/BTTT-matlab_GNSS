function [icol] = PickEPH_GLO(eph_glo,prn,time)

% function [icol] = PickEPH_GLO2(eph_glo,prn,time)
% 
%  DO: GLONASS ��۱˵��� ��Ŀ��� �Էµ� �ð��� ���� ���� ���������� ��.
% 
%  - Coded by Mi-So Kim, June 27th 2014 


icol = 0;
isat = find(eph_glo(:,1) == prn);

n = length(isat);

if n == 0
    return
end;

icol = isat(1);
dtmin = eph_glo(icol,2) - time;

for k = 2:n
   kk = isat(k);
   dt = eph_glo(kk,2) - time;
   if dt < 0
       if abs(dt) < abs(dtmin)
       icol = kk;
       dtmin = dt;
       end
   end   
end