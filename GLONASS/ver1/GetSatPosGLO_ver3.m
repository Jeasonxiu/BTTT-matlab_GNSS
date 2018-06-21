function [sat_pos,sat_vel,sat_LS] = GetSatPosGLO_ver3(SatPosArr_before,jcol,tc,deltat)
%
% <input>
%   eph_glo = GLONASS BRDC array
%      icol = PickEPH_GLO�� output
%       tc  = ��ȣ���� �ӵ��� ����� �ð� (ts = gs - TTs)
%    deltat = ���� ����
%
% <output>
%   sat_pos = �ش� �ð������� GLONASS�� ��ġ
%   sat_vel = �ش� �ð������� GLONASS�� �ӵ�
%
%  Copyright: Mi-so Kim, January 14th, 2015
%
%%
gs = SatPosArr_before(jcol,2); x = SatPosArr_before(jcol,1:11);
%%
% WGS_xyz = PZ2WGS(x(3:5)); % ��ǥ��ȯ - ��ǥ��ȯ�� ��� �κп� ���������� ���/ 1)��ǥ��� �� ��ǥ��ȯ 2) ��ǥ��ȯ �� ��ǥ���
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
%         for k = 1:10  % �񱳰��, �ݺ������ ��������� �̹��Ͽ� �������� ����
    [xyz,vxyz,xyzLS] = RK4_my(x(a,:), hat);
%         end
    sat_pos = xyz;
    sat_vel = vxyz;
    sat_LS = xyzLS;
end