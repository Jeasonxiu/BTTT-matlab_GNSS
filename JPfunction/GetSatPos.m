function [SatPos] = GetSatPos(eph, i_eph, t_GPS)
%
%function [SatPos] = GetSatPos(eph, i_eph, t_GPS)
%
%   Computes the XYZ coordinates of a given GPS satellite
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
i_dot     = eph(i_eph,16);
week_no   = eph(i_eph,17);

mu = 3.986005e14;
Omega_dot_e = 7.2921151467e-5;

A = sqtA.^2;
n_0 = sqrt(mu./(A.^3));
n = n_0 + dn;

t_oe_f = t_oe;

% t_off = 0.075;

%== t_rec is the time of signal reception referenced to time of ephemerids
t_rec = t_GPS - t_oe_f;

%== set t_rec to t_k
t_k = t_rec;

%== Mean Motion
M_k = M_0 + n.*t_k;

%== translate mean anomaly into eccentric anomaly
E_k = ecce_anom(M_k,e,3);

%== translate eccentric anomaly into true anomaly
nu_k = atan((sqrt(1-e.^2).*sin(E_k))./(cos(E_k)-e));
%== atan gives angles only  between -pi/2 and pi/2, compensate for this 
ind_nu_k = find(cos(E_k)-e<0);
if ~isempty(ind_nu_k)
  nu_k(ind_nu_k) = nu_k(ind_nu_k) + pi;
end

%== Argument of latitude
Phi_k = nu_k + omega;

du_k = C_us.*sin(2*Phi_k) + C_uc.*cos(2*Phi_k);  %Second harmonic perturbations
dr_k = C_rs.*sin(2*Phi_k) + C_rc.*cos(2*Phi_k);
di_k = C_is.*sin(2*Phi_k) + C_ic.*cos(2*Phi_k);

u_k = Phi_k + du_k;                % Corrected arg. of lat.
r_k = A.*(1-e.*cos(E_k)) + dr_k;   % Corrected radius
i_k = i_0+di_k + i_dot.*t_k;       % Corrected inclination

%== Sat pos. in orb. plane
xp_k = r_k.*cos(u_k);
yp_k = r_k.*sin(u_k);

%== The corrected long. of asc. node in the inertial reference frame
Omega_k = Omega_0 + (Omega_dot-Omega_dot_e).*t_rec - Omega_dot_e.*t_oe;
 
%== Sat pos. in the inertial coord system coinsiding with ITRF at t_rec
x_k = xp_k.*cos(Omega_k) - yp_k.*cos(i_k).*sin(Omega_k);
y_k = xp_k.*sin(Omega_k) + yp_k.*cos(i_k).*cos(Omega_k);
z_k = yp_k.*sin(i_k);

p_sat=[x_k y_k z_k];

SatPos = p_sat';
