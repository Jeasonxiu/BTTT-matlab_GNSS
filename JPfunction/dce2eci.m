function C_ecef2eci = dce2eci(rcvrTimeWeek, rcvrTimeSec, rcvrTime2gpsTime)
%
% Generate DCM from ECEF (Earth-Centered Earth-Fixed) to ECI (Earth-Centered Inertial) coordinate systems
%   C_ecef2eci = dce2eci(rcvrTimeWeek, rcvrTimeSec, rcvrTime2gpsTime)
% 
% Inputs:
%   rcvrTimeWeek     - number of receiver weeks since GPS epoch
%   rcvrTimeSec      - number of receiver seconds into GPS week
%   rcvrTime2gpsTime - conversion from receiver time to GPS time
%
% Output:
%   C_ecef2eci - DCM to convert a vector from ECEF coordinates to ECI coordinates
%
% Example:   v_eci = C_ecef2eci * v_ecef, 
%                    C_ecef2eci' * v_eci = v_ecef
%
% Author(s): Steve Young, May 2006
% Copyright 2006, NAVSYS Corporation

deg2rad = pi/180;  % conversion from degrees to radians

% Specify timezone as integer hour offset to UTC (we want UTC)
timezone = 0;

% Compute UTC time from GPS-time and week number
gpsTime = rcvrTimeSec + rcvrTime2gpsTime;
UTC = fromgpst(rcvrTimeWeek, gpsTime, timezone, [], 0); % [year,month,day,hours,minutes,seconds]

% Compute Modified Julian Date (# days since 01 January 1972 00:00:00)
MJD = julian(UTC(1), UTC(2), UTC(3));
MJD = MJD + (UTC(4)/24) + (UTC(5)/1440) + (UTC(6)/86400);

% Compute number of Julian Centuries since J2000 epoch (01 January 2000 12:00:00)
%   Reference: equation 3-40 on pg. 188 of Vallado
T_UT1 = (MJD - 51544.5) / 36525;

T2 = T_UT1^2;
T3 = T_UT1^3;
T4 = T_UT1^4;

% Compute rotation angles representing the combined effects of general
% precession, to rotate between the ECI (FK5 J2000) frame and the Mean Equinox
% of Date (MOD) frame (from Vallado, pg. 215 equations 3-56 and 3-57):
%       [1] zeta  = rotation along the mean equator at J2000
%       [2] theta = rotation to the mean equator of date
%       [3] z     = rotation to the equinox of date
%
% The rotation for transformations from ECI (FK5 J2000) to MOD is given by:
%       v_mod = ROT3(-z)*ROT2(theta)*ROT3(-zeta)*v_eci

zeta  = (0.6406161*T_UT1 + 0.0000839*T2 +  5.0e-6*T3)*deg2rad;
theta = (0.5567530*T_UT1 - 0.0001185*T2 - 1.16e-5*T3)*deg2rad;
z     = (0.6406161*T_UT1 + 0.0003041*T2 +  5.1e-6*T3)*deg2rad;

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
% frame and the True of Date (TOD) frame (from Vallado, pg. 217-219):
%       [1] dPsi = nutation in longitude (from mean equator of date)
%       [2] dEps = nutation in obliquity (from mean equator of date)
%           eps  = true obliquity of the ecliptic
%           epsBar = mean obliquity of the ecliptic
%
% The rotation for transformations from MOD to TOD is given by:
%       v_tod = ROT1(-eps)ROT3(-dPsi)ROT1(epsBar)*v_mod

epsBar = (23.439291 - 0.0130042*T_UT1 - 1.64e-7*T2 + 5.04e-7*T3)*deg2rad;

% Now compute the 5 values needed in the dPsi & dEps equations
% Mean Anomaly for the Moon:
M_moon  = (134.96340251 + 198.8675605*T_UT1 + 0.0088553*T2 + 1.4343e-5*T3 ...
           - 6.797e-6*T4 + 1325*360*T_UT1)*deg2rad;

% Mean Anomaly for the Sun:
M_sun   = (357.52910918 + 359.0502911*T_UT1 - 0.0001537*T2 - 3.8e-8*T3 ...
           - 3.190e-9*T4 + 99*360*T_UT1)*deg2rad;

% Mean Argument of Latitude of the Moon
uM_moon = ( 93.27209062 +  82.0174577*T_UT1 - 0.0035420*T2 + 2.88e-7*T3 ...
           + 1.160e-9*T4 + 1342*360*T_UT1)*deg2rad;
       
% Mean Elongation from the Sun
D_sun   = (297.85019547 + 307.1114469*T_UT1 - 0.0017696*T2 + 1.831e-6*T3 ...
           - 8.800e-9*T4 + 1236*360*T_UT1)*deg2rad;
       
% Longitude of the Ascending Node of the Mean Lunar Orbit:
Omega_moon = (125.04455501 - 134.1361851*T_UT1 + 0.0020756*T2 + 2.139e-6*T3 ...
              - 1.650e-8*T4 - 5*360*T_UT1)*deg2rad;

% Get the coefficients for the trigonometric series approximation
[Ai, Bi, Ci, Di, ai] = nutation_coef;

% Compute the series using the first 20 terms (12 May 2006)
%    dPsi = SUM{ (Ai + Bi*T_tdb)*sin(ap) }
%    dEps = SUM{ (Ci + Di*T_tdb)*cos(ap) }
%    ap = a1i*M_moon + a2i*M_sun + a3i*uM_moon + a4i*D_sun + a5i*Omega_moon

ap   = ai * [ M_moon; M_sun; uM_moon; D_sun; Omega_moon ];
dPsi = sum( (Ai + Bi*T_UT1) .* sin(ap) );
dEps = sum( (Ci + Di*T_UT1) .* cos(ap) );

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
% and the Pseudo-Earth-Fixed (PEF) frame (from Vallado, pg. 219):
%       [1] theta_GMST = Greenwich Mean Sidereal Time
%       [2] EQ_equinox = Equation of the Equinoxes; difference between the
%                        mean and true equinoxes projected onto the true equator
%       [3] theta_AST  = Greenwich Apparent Sidereal Time
%
%       where:  theta_AST = theta_GMST + EQ_equinox;
%
% The rotation for transformations from TOD to PEF is given by:
%       v_pef = ROT3(theta_AST)*v_tod

% Determine Greenwich Mean Sidereal Time in seconds and convert to radians (within 2*pi of zero)
%   Reference: equation 3-45 on pg. 191 of Vallado
theta_GMST_seconds = 67310.54841 + (3155760000 + 8640184.812866)*T_UT1 + 0.093104*T2 - 6.2e-6*T3;
theta_GMST = rem(theta_GMST_seconds, 86400)*(pi/43200);  % in radians

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
% and the ECEF (WGS84) frame (from Vallado, pg. 220):
%       [1] xp = displacement of CEP pole relative to IRP pole, measured
%                south along the 0 deg longitude meridian
%       [2] yp = displacement of CEP pole relative to IRP pole, measured
%                south along the 90 deg west longitude meridian
%
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