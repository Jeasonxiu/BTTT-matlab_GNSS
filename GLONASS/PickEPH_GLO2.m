function [icol] = PickEPH_GLO2(eph_glo,prn,time)
% eph_glo = ephGLO;
% prn = 15;
% time = 5220;
% 2014.06.27 미소
% eph_glo array에서 입력된 time에 따른 행을 pick.


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