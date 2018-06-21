function [dRel] = GetRelBRDC(eph, i_eph, gs)
%
% DO: Extract the relativistic effect form BRDC information
%
% <input>   eph: Ephermeris array obtained from ReadEPH
%           i_eph: Index for ephemeris
%           gs: Epoch time to get the relativistic effect
%
% <output>  dRel: Relativity effect
%
% Copyright: Kwan-Dong Park, December 31, 2013 @LDEO
%
% �μ� �ڸ�Ʈ: ��뼺ȿ�� ��꿡�� e, sqrt(semi-major axis), Eccentric anomaly�� �ʿ���

%% EPH���� ������ ������ - ���ʿ��� �͵��� ��������
e         = eph(i_eph,05);
sqtA      = eph(i_eph,07);
dn        = eph(i_eph,02);
t_oe      = eph(i_eph,08);
M_0       = eph(i_eph,03);
%% ����� ��� mu�� �ʿ���
mu = 3.986005e14;
%% ��� ����
A = sqtA^2;
n_0 = sqrt(mu/(A^3));
n = n_0 + dn;
t_k = gs - t_oe;
M_k = M_0 + n*t_k;
E_k = ecce_anom(M_k, e, 20);
%% ��뼺 ȿ�� ���
result=[A,n_0,n,M_k,E_k];
dRel = -4.442807633*10^(-10)*e*sqtA*sin(E_k);
