% function [r_eci,v_eci] = ecef2eci(jd,r_ecef,v_ecef)
%
% function [r,v] = ecef2eci(jd,r_ecef,v_ecef)
%
%   Determining Keplerian Elements from Satellites's POS, VEL
%   vector(ecef:WGS84)
%   
%   input jd[1 x 1] : Sidereal time to theta
%   input r_ecef[1 x 3] : Position Vector(ECEF:WGS84)
%   input v_ecef[1 x 3] : Velocity vector(ECEF:WGS84)
%
%   output r_eci[1 x 3] : Position Vector(ECI)
%   output v_eci[1 x 3] : Velocity vector(ECI)
%
%   Example : [r,v] = ecef2eci(jd,r_ecef,v_ecef)
%
%   coded by Joonseong Gim, April 8, 2016
%

%% gw, gs to Julian day
clear all; close all;
gw = 1803;
gs = 346500;
[jd] = gwgs2jd(gw, gs);
load('glomat.mat');
JD = EphGlo(1,18);
R = EphGlo(1,3:5);                                              % Position Vector(PZ-90)
V = EphGlo(1,6:8);                                              % Velocity Vector(PZ-90)
r_ecef = PZ2WGS(R);                                             % Position Vector(WGS84)
v_ecef = PZ2WGS(V);                                             % Velocity Vector(WGS84)

% conversion from degrees to radians
deg2rad = pi/180;
% Compute number of Julian Centuries since J2000 epoch (01 January 2000 12:00:00)
Ttdb = (jd -2451545)/36525;
T2 = Ttdb^2; 
T3 = Ttdb^3;
T4 = Ttdb^4;

% Compute rotation angles representing the combined effects of general
% precession, to rotate between the ECI (FK5 J2000) frame and the Mean Equinox
% of Date (MOD) frame (from Vallado, pg. 215 equations 3-56 and 3-57):
zeta  = (0.6406161*Ttdb + 0.0000839*T2 +  5.0e-6*T3)*deg2rad;
theta = (0.5567530*Ttdb - 0.0001185*T2 - 1.16e-5*T3)*deg2rad;
z     = (0.6406161*Ttdb + 0.0003041*T2 +  5.1e-6*T3)*deg2rad;

% The rotation for transformations from ECI (FK5 J2000) to MOD is given by:
%       v_mod = ROT3(-z)*ROT2(theta)*ROT3(-zeta)*v_eci
cZeta  = cos(zeta);
sZeta  = sin(zeta);
cTheta = cos(theta);
sTheta = sin(theta);
cZ     = cos(z);
sZ     = sin(z);

C_eci2mod = [ cTheta*cZ*cZeta - sZ*sZeta,  -sZeta*cTheta*cZ - sZ*cZeta,  -sTheta*cZ; ...
              sZ*cTheta*cZeta + sZeta*cZ,  -sZ*sZeta*cTheta + cZ*cZeta,  -sTheta*sZ; ...
                            sTheta*cZeta,                -sTheta*sZeta,      cTheta ];
                        
% Compute rotation angles representing the periodic effects of nutation
% (primarily from the Moon), to rotate between the Mean of Date (MOD)
epsBar = (23.439291 - 0.0130042*Ttdb - 1.64e-7*T2 + 5.04e-7*T3)*deg2rad;

% Now compute the 5 values needed in the dPsi & dEps equations
% Mean Anomaly for the Moon:
M_moon  = rem((134.96340251 + 198.8675605*Ttdb + 0.0088553*T2 + 1.4343e-5*T3 ...
           - 6.797e-6*T4 + 1325*360*Ttdb),360)*deg2rad;

% Mean Anomaly for the Sun:
M_sun   = rem((357.52910918 + 359.0502911*Ttdb - 0.0001537*T2 - 3.8e-8*T3 ...
           - 3.190e-9*T4 + 99*360*Ttdb),360)*deg2rad;

% Mean Argument of Latitude of the Moon
uM_moon = rem((93.27209062 +  82.0174577*Ttdb - 0.0035420*T2 + 2.88e-7*T3 ...
           + 1.160e-9*T4 + 1342*360*Ttdb),360)*deg2rad;
       
% Mean Elongation from the Sun
D_sun   = rem((297.85019547 + 307.1114469*Ttdb - 0.0017696*T2 + 1.831e-6*T3 ...
           - 8.800e-9*T4 + 1236*360*Ttdb),360)*deg2rad;
       
% Longitude of the Ascending Node of the Mean Lunar Orbit:
Omega_moon = rem((125.04455501 - 134.1361851*Ttdb + 0.0020756*T2 + 2.139e-6*T3 ...
              - 1.650e-8*T4 - 5*360*Ttdb),360)*deg2rad;
          
% Get the coefficients for the trigonometric series approximation
[Ai, Bi, Ci, Di, ai] = nutation_coef;

% Compute the series using the first 20 terms (12 May 2006)
%    dPsi = SUM{ (Ai + Bi*T_tdb)*sin(ap) }
%    dEps = SUM{ (Ci + Di*T_tdb)*cos(ap) }
%    ap = a1i*M_moon + a2i*M_sun + a3i*uM_moon + a4i*D_sun + a5i*Omega_moon

ap   = ai * [ M_moon; M_sun; uM_moon; D_sun; Omega_moon ];
dPsi = sum( (Ai*0.0001 + Bi*0.0001*Ttdb) .* sin(ap) )/3600;
dEps = sum( (Ci*0.0001 + Di*0.0001*Ttdb) .* cos(ap) )/3600;

eps = epsBar + dEps;

% Define the DCM from MOD to TOD using the angles computed above
c_dPsi = cos(dPsi);
s_dPsi = sin(dPsi);
c_eBar = cos(epsBar);
s_eBar = sin(epsBar);
c_eps  = cos(eps);
s_eps  = sin(eps);

C_mod2tod = [       c_dPsi,                      -s_dPsi*c_eBar,                      -s_dPsi*s_eBar; ...
              s_dPsi*c_eps,  c_eps*c_dPsi*c_eBar + s_eps*s_eBar,  s_eBar*c_eps*c_dPsi - s_eps*c_eBar; ...
              s_eps*s_dPsi,  s_eps*c_dPsi*c_eBar - s_eBar*c_eps,  s_eps*s_eBar*c_dPsi + c_eps*c_eBar ];
          
% Compute rotation angles representing Greenwich apparent sidereal time
% (due to Earth's rotation), to rotate between the True of Date (TOD) frame
% The rotation for transformations from TOD to PEF is given by:
%       v_pef = ROT3(theta_AST)*v_tod

% Determine Greenwich Mean Sidereal Time in seconds and convert to radians (within 2*pi of zero)
%   Reference: equation 3-45 on pg. 191 of Vallado
theta_GMST_seconds = 67310.54841 + (3155760000 + 8640184.812866)*Ttdb + 0.093104*T2 - 6.2e-6*T3;
theta_GMST = rem(theta_GMST_seconds, 86400)/240;
% theta_GMST_seconds = (100.4606184+36000.77005361*Ttdb + 0.00038793*T2 - 2.6e-8*T3)/3600;
% theta_GMST = rem(theta_GMST_seconds, 24)*deg2rad;  % in radians

% Compute the Equation of the Equinoxes (in radians)
%   Reference: equation 3-65 on pg. 219 of Vallado
EQ_equinox = dPsi*cos(eps) + 1.28e-8*sin(Omega_moon) + 3.054e-10*sin(2*Omega_moon);

% Finally, compute Greenwich Apparent Sidereal Time
theta_AST = theta_GMST + EQ_equinox;
c_thetaAST = cos(theta_AST);
s_thetaAST = sin(theta_AST);

% Construct the DCM for rotation from TOD to PEF
C_tod2pef = [ c_thetaAST, s_thetaAST,   0; ...
             -s_thetaAST, c_thetaAST,   0; ...
                       0,          0,   1 ];
                   
% Compute rotation angles representing polar motion (with respect to the 
% Earth's crust), to rotate between the Pseudo-Earth-Fixed (PEF) frame 

% The rotation for transformations from PEF to ECEF is given by:
%       v_ecef = ROT2(-xp)ROT1(-yp)*v_pef
                   
xp = 0;  % Neglect polar motion for now (5/15/06)
yp = 0;  % Neglect polar motion for now (5/15/06)

% Use the small-angle approximation DCM for rotation from PEF to ECEF
C_pef2ecef = [  1,  0,  xp; ...
                0,  1, -yp; ...
              -xp, yp,   1 ];
          
% Combine the rotations computed above to generate a DCM rotation from ECEF
% to ECI:  v_eci = C_ecef2eci*v_ecef = {C_pef2ecef*C_tod2pef*C_mod2tod*C_eci2mod}'*v_ecef
C_eci2ecef = C_pef2ecef * C_tod2pef * C_mod2tod * C_eci2mod;
C_ecef2eci = C_eci2ecef';
r_eci = C_ecef2eci*r_ecef'

