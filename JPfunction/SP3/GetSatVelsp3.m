function [vel] = GetSatVelsp3(sp3, prn, tc)
%
%function [vel] = GetSatVelsp3(sp3, prn, tc)
%
% DO: Compute the satellite velocity based on SP3
%
% <input>   sp3: SP3 array
%           prn: PRN number 
%           tc: time-epoch corrected for travel time
%
% <output>  vel: Velocity Vector [3x1 row vector]
%
% Copyright: Kwan-Dong Park, January 21, 2014 @LDEO
%

%% �Լ� �ۼ� �� �Է� �׽�Ʈ
% prn = 1;
% tc = 15.5;
%% tc1�� tc2���� pos1/pos2�� ����ϰ� �̸� �������� vel(�ӵ�) ���
dt = 1e-3;
tc1 = tc;
tc2 = tc + dt;
out1 = IntpSP3e1(sp3, prn, tc1); pos1 = out1(2:4);
out2 = IntpSP3e1(sp3, prn, tc2); pos2 = out2(2:4);
vel = (pos2 - pos1)/dt;
%% ����� ���� ��� Ȯ��
% fprintf('%13.3f %13.3f %13.3f\n', pos1(1:3))
% fprintf('%13.3f %13.3f %13.3f\n', pos2(1:3))
% fprintf('%13.3f %13.3f %13.3f %13.3f\n', vel(1:3), norm(vel))