function [icol] = PickEPH_GLO2(eph_glo,prn,time)

% function [icol] = PickEPH_GLO2(eph_glo,prn,time)
% 
%  DO: GLONASS 방송궤도력 행렬에서 입력된 시간에 가까운 행을 가져오도록 함.
% 
%  - Coded by Mi-So Kim, June 27th 2015 


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
%    if dt < 0
       if abs(dt) <= abs(dtmin)
       icol = kk;
       dtmin = dt;
       end
%    end   
end