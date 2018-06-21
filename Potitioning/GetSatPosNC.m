function [SatPos] = GetSatPosNC(eph, i_eph, gs)
%
%function [SatPos] = GetSatPosNC(eph, i_eph, gs)
%
%** 신호전달시간을 보정하지 않은 gs(gps week second)를 그대로 사용하고, 
%   신호전달시간동안의 지구자전도 고려하지 않고 위성좌표를 계산한 다음
%   SP3의 참값과 비교하기 위한 코드임.
%

C_rs      = eph(i_eph,01);
dn        = eph(i_eph,02);
M_0       = eph(i_eph,03);
C_uc      = eph(i_eph,04);
e         = eph(i_eph,05);
C_us      = eph(i_eph,06);
sqtA      = eph(i_eph,07);
t_oe      = eph(i_eph,08);
C_ic      = eph(i_eph,09);
Omega_0   = eph(i_eph,10);
C_is      = eph(i_eph,11);
i_0       = eph(i_eph,12);
C_rc      = eph(i_eph,13);
omega     = eph(i_eph,14);
Omega_dot = eph(i_eph,15);
i_dot     = eph(i_eph,16);      % week_no   = eph(i_eph,17); 현재로선 불필요한 부분
%% WGS84에서 정의된 상수
mu = 3.986005e14;
Omega_dot_e = 7.2921151467e-5;
%% 이후 코드에서 사용하게 될 t_k 정의: 단순히 gs에서 t_eo를 뺌
t_k = gs - t_oe;
%% 위성궤도 연산 시작(참고문헌은 다양함, 여기서 굳이 코드 설명을 붙일 필요가 없다고 봄 
A = sqtA^2;
n_0 = sqrt(mu / (A^3));
n = n_0 + dn;
M_k = M_0 + n*t_k;
E_k = ecce_anom(M_k, e, 10);
nu_k = atan2((sqrt(1 - e^2)*sin(E_k)), (cos(E_k) - e));
Phi_k = nu_k + omega;
du_k = C_us*sin(2*Phi_k) + C_uc*cos(2*Phi_k);  
dr_k = C_rs*sin(2*Phi_k) + C_rc*cos(2*Phi_k);
di_k = C_is*sin(2*Phi_k) + C_ic*cos(2*Phi_k);
u_k = Phi_k + du_k;                 
r_k = A*(1 - e*cos(E_k)) + dr_k;    
i_k = i_0 + di_k + i_dot*t_k;       
xp_k = r_k*cos(u_k);                
yp_k = r_k*sin(u_k);
Omega_k = Omega_0 + (Omega_dot - Omega_dot_e)*t_k - Omega_dot_e*t_oe;
x_k = xp_k*cos(Omega_k) - yp_k*cos(i_k)*sin(Omega_k);
y_k = xp_k*sin(Omega_k) + yp_k*cos(i_k)*cos(Omega_k);
z_k = yp_k*sin(i_k);
SatPos = [x_k y_k z_k];