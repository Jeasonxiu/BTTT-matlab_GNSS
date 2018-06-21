function [f] = EOM_glo(x)
% 
% function [f] = EOM_glo(x)
% 
% Copyright: Kwan-Dong Park & Mi-So Kim , January 14th, 2015
%
%--- Modifications --
% 1/14/2015: 교수님께서 작성하신 기본틀을 행렬순서에 맞게 수정함

%% 상수 정의 --- 재검토 필요(ICD에서 제시한 상수랑 다른 수치가 있음) 1/2/2015; GLONASS ICD 2008 확인
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
%% 반경(radius) 계산
r = sqrt(x(3)^2 + x(4)^2 + x(5)^2);
r2 = r*r;
r3 = r2*r;
r5 = r3*r2;
%% 운동방정식 - 반복되는 항을 미리 정의해서 연산 속도 높임
mur3 = mu/r3;
c20p = (3/2)*c20*mu*aE2/r5;
z2r2 = x(5)^2/r2; % Z^2/r2
%% 운동방정식 - GLONASS 운동 방정식 ---?? Luni-Solar는 어디서 구하나? 1/2/2014 - ans) GLONASS 방송궤도력에서 제공
LS_X = x(9); LS_Y = x(10); LS_Z = x(11);
f(1) = -mur3*x(3) + c20p*x(3)*(1 - 5*z2r2) + LS_X + om2*x(3) + 2*om*x(7);   % (dVx/dt): K(X)
f(2) = -mur3*x(4) + c20p*x(4)*(1 - 5*z2r2) + LS_Y + om2*x(4) - 2*om*x(6);   % (dVy/dt): K(Y)
f(3) = -mur3*x(5) + c20p*x(5)*(3 - 5*z2r2) + LS_Z;                          % (dVz/dt): K(Z)