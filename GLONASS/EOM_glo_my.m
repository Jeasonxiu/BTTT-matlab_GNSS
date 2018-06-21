function [f] = EOM_glo_my(x,acc)
% 
% function [f] = EOM_glo(x)
% <input>:
%  x = initial value          ex) [3x 4y 5z 6xdot 7ydot 8zdot] 
%  acc = initial value        ex) [ 1xLS 2yLS 3zLS ]
% 
% <output>: glonass orbit differential equations
% f [vel_x vel_y vel_z dx dy dz]
% Copyright: Kwan-Dong Park & Mi-So Kim , January 14th, 2015
%            Jong-Seok Kim , July 7, 2016
%--- Modifications --
% 1/14/2015: �����Բ��� �ۼ��Ͻ� �⺻Ʋ�� ��ļ����� �°� ������
% 7/01/2016: input ����

%% ��� ���� --- ����� �ʿ�(ICD���� ������ ����� �ٸ� ��ġ�� ����) 1/2/2015; GLONASS ICD 2008 Ȯ��
mu = 3.9860044109e14;        %: mu = GM
c20 = -sqrt(5)*484165e-9;   %: c20 = -J2
aE = 6378136;               %: aE = Mean Earth radius
aE2 = aE*aE;
om = 7.2921151467e-5;       %: om = omega
om2 = om*om;
% mu = 398600.4418e9;        %: mu = GM
% c20 = -1082625.75e-9;   %: c20 = -J2
% aE = 6378136;               %: aE = Mean Earth radius
% aE2 = aE*aE;
% om = 7.292115e-5;       %: om = omega
% om2 = om*om;
%% �ݰ�(radius) ���
r = sqrt(x(1)^2 + x(2)^2 + x(3)^2);
r2 = r*r;
r3 = r2*r;
r5 = r3*r2;
%% ������� - �ݺ��Ǵ� ���� �̸� �����ؼ� ���� �ӵ� ����
mur3 = mu/r3;
c20p = (3/2)*c20*mu*aE2/r5;
z2r2 = x(3)^2/r2; % Z^2/r2
%% ������� - GLONASS � ������ ---?? Luni-Solar�� ��� ���ϳ�? 1/2/2014 - ans) GLONASS ��۱˵��¿��� ����
LS_X = acc(1); LS_Y = acc(2); LS_Z = acc(3);
f(1) = x(4);
f(2) = x(5);
f(3) = x(6);
f(4) = -mur3*x(1) + c20p*x(1)*(1 - 5*z2r2) + LS_X + om2*x(1) + 2*om*x(5);   % (dVx/dt): K(X)
f(5) = -mur3*x(2) + c20p*x(2)*(1 - 5*z2r2) + LS_Y + om2*x(2) - 2*om*x(4);   % (dVy/dt): K(Y)
f(6) = -mur3*x(3) + c20p*x(3)*(3 - 5*z2r2) + LS_Z;                          % (dVz/dt): K(Z)
