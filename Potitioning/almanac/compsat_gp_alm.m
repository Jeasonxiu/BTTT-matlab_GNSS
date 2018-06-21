function [p_sat] = compsat_gp_alm(t_GPS, alm, i_alm)

%
%function [p_sat] = compsat_gp_alm(t_GPS, alm, i_alm)
%
% Computes the XYZ coordinates of a given GPS satellite
%
% July 04, 2008
%
% [Hye-In Kim, Inha University]

PRN = alm(i_alm,01);
Health = alm(i_alm,02);
e = alm(i_alm,03);
t_oa = alm(i_alm,04);
i = alm(i_alm,05);
Omega_dot = alm(i_alm,06);
sqtA = alm(i_alm,07);
Omega_0 = alm(i_alm,08);
omega = alm(i_alm,09);
M_0 = alm(i_alm,10);
Af0 = alm(i_alm,11);
Af1 = alm(i_alm,12);
week = alm(i_alm,13);

if M_0 > 3.1416
    M_0 = M_0 - 6.2832
end

Omega_dot_e = 7.2921151467 * 10^(-5);
mu = 3.986004418 * 10^14;
a = sqtA.^2;
n_0 = sqrt(mu./(a.^3)); 


% t_k = t_GPS - (t_oa-(86400*2));
t_k = t_GPS - t_oa
% if t_k > 302400
%     t_k = t_k - 604800;
% elseif t_k < -302400
%     t_k = t_k + 604800;
% end

M_k = M_0 + (n_0 * t_k);
E_k = Kepler_solve(M_k,e,5);
nu_k = atan2((sqrt(1-e.^2)*sin(E_k)),(cos(E_k)-e)); 

p = a*(1-e^2);
r = p/(1+e*cos(nu_k));
x = r*cos(nu_k);
y = r*sin(nu_k);

Omega_k = Omega_0+((Omega_dot-Omega_dot_e).*t_k)-(Omega_dot_e.*t_oa);

T1 = [cos(Omega_k) sin(Omega_k) 0; -sin(Omega_k) cos(Omega_k) 0; 0 0 1];
T2 = [1 0 0;0 cos(i) sin(i); 0 -sin(i) cos(i)];
T3 = [cos(omega) sin(omega) 0; -sin(omega) cos(omega) 0; 0 0 1];

T = T3 * T2 * T1;

[p_sat] = T'*[x; y; 0];
p_sat = p_sat';
   

