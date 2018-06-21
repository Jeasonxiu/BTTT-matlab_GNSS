function [SatPos] = GetSatPos_GC_almanac0(alm, prn, gs)
%
%function [SatPos] = GetSatPos_GC_almanac(alm, prn, gs)
%
%** �˸��� �����͸� �̿��� ������ġ ����
%
%

% clear all
% close all
%
% load('almanac_17255.mat');
%
% alm = [alm_gps;alm_bds];
% prn = 101;
% gs = 440000;

e         = alm(find(alm(:,1) == prn),3);
t_oa      = alm(find(alm(:,1) == prn),4);
i_0       = alm(find(alm(:,1) == prn),5);
Omega_dot = alm(find(alm(:,1) == prn),6);
sqtA      = alm(find(alm(:,1) == prn),7);
Omega_0   = alm(find(alm(:,1) == prn),8);
omega     = alm(find(alm(:,1) == prn),9);
M_0       = alm(find(alm(:,1) == prn),10);
Af0       = alm(find(alm(:,1) == prn),11);
Af1       = alm(find(alm(:,1) == prn),12);
week_no   = alm(find(alm(:,1) == prn),13);      % week_no   = eph(i_eph,17); ����μ� ���ʿ��� �κ�


%% WGS84���� ���ǵ� ���
mu = 3.986005e14;
Omega_dot_e = 7.2921151467e-5;
%% ���� �ڵ忡�� ����ϰ� �� t_k ����: �ܼ��� gs���� t_oa�� ��
t_k = gs - t_oa;
%% �����˵� ���� ����(�������� �پ���, ���⼭ ���� �ڵ� ������ ���� �ʿ䰡 ���ٰ� ��
A = sqtA^2;
n_0 = sqrt(mu / (A^3));
M_k = M_0 + n_0*t_k;
E_k = ecce_anom(M_k, e, 3);
nu_k = atan2((sqrt(1 - e^2)*sin(E_k)), (cos(E_k) - e));

p = A*(1-e^2);
r = p/(1+e*cos(nu_k));
x = r*cos(nu_k);
y = r*sin(nu_k);

Omega_k = Omega_0+((Omega_dot-Omega_dot_e).*t_k)-(Omega_dot_e.*t_oa);
if prn < 200
    T1 = [cos(Omega_k) sin(Omega_k) 0; -sin(Omega_k) cos(Omega_k) 0; 0 0 1];
    T2 = [1 0 0;0 cos(i_0) sin(i_0); 0 -sin(i_0) cos(i_0)];
    T3 = [cos(omega) sin(omega) 0; -sin(omega) cos(omega) 0; 0 0 1];
elseif prn > 200
    if prn > 205
        i = 0.3*pi + i_0;
    else
        i = 0.00*pi + i_0;
    end
    T1 = [cos(Omega_k) sin(Omega_k) 0; -sin(Omega_k) cos(Omega_k) 0; 0 0 1];
    T2 = [1 0 0;0 cos(i) sin(i); 0 -sin(i) cos(i)];
    T3 = [cos(omega) sin(omega) 0; -sin(omega) cos(omega) 0; 0 0 1];
end
T = T3 * T2 * T1;


[SatPos] = T'*[x; y; 0];
SatPos = SatPos';
