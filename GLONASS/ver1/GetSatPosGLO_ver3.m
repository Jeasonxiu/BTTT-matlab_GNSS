function [sat_pos,sat_vel,sat_LS] = GetSatPosGLO_ver3(SatPosArr_before,jcol,tc,deltat)
%
% <input>
%   eph_glo = GLONASS BRDC array
%      icol = PickEPH_GLO의 output
%       tc  = 신호전달 속도가 적용된 시간 (ts = gs - TTs)
%    deltat = 적분 간격
%
% <output>
%   sat_pos = 해당 시각에서의 GLONASS의 위치
%   sat_vel = 해당 시각에서의 GLONASS의 속도
%
%  Copyright: Mi-so Kim, January 14th, 2015
%
%%
gs = SatPosArr_before(jcol,2); x = SatPosArr_before(jcol,1:11);
%%
% WGS_xyz = PZ2WGS(x(3:5)); % 좌표변환 - 좌표변환을 어디 부분에 적용할지는 고민/ 1)좌표계산 후 좌표변환 2) 좌표변환 후 좌표계산
% x(3:5) = WGS_xyz;
%%
dif = tc - gs;
if dif < 0 % backward
    deltat = -deltat;
    it = fix(dif/deltat);
    hat = rem(dif,deltat);
else % forward
    it = fix(dif/deltat);
    hat = rem(dif,deltat);
end
%%
a = 1;
if it == 0
%         for k = 1:10
    [xyz,vxyz,xyzLS] = RK4_my(x, hat);
%         end
    sat_pos = xyz;
    sat_vel = vxyz;
    sat_LS = xyzLS;
    return;
else
    it = abs(it);
    for i = 1:it
        [xyz,vxyz] = RK4_my(x(i,:), deltat);
        a = a + 1;
        x(a,2) = x(a-1,2) + deltat;
        x(a,3:5) = xyz;
        x(a,6:8) = vxyz;
        x(a,9:11) = x(a-1,9:11);
    end
%         for k = 1:10  % 비교결과, 반복계산의 향상정도가 미미하여 적용하지 않음
    [xyz,vxyz,xyzLS] = RK4_my(x(a,:), hat);
%         end
    sat_pos = xyz;
    sat_vel = vxyz;
    sat_LS = xyzLS;
end