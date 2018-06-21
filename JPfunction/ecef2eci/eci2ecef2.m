function [recef,vecef] = eci2ecef2(jd,r_eci,v_eci)

% function [r_ecef,v_ecef] = eci2ecef2(jd,r_eci,v_eci)
%
%   Determining Keplerian Elements from Satellites's POS, VEL
%   vector(ecef:WGS84)
%   
%   input jd[1 x 1] : Julian day include(dAT, dUT1)
%   input r_eci[1 x 3] : Position Vector(ECI)
%   input v_eci[1 x 3] : Velocity vector(ECI)
%
%   output r_ecef[3 x 1] : Position Vector(ECEF:WGS84)
%   output v_ecef[3 x 1] : Velocity vector(ECEF:WGS84)
%
%   Example : [r_ecef,v_ecef] = eci2ecef2(jd,r_eci,v_eci)
%
%   coded by Joonseong Gim, April 15, 2016
%
%

%% 테스트를 위한 파라미터 입력
% clear all; close all;
% r_ecef = [-9529416.99219000,-15029482.9102000,-18259631.8359000];              % example : p.230, vallado 4th edition
% v_ecef = [-198.983192444000,-2495.43666840000,2157.84358978000];               % example : p.230, vallado 4th edition
% jd=[2456869.51041667];                                           % example : p.230, vallado 4th edition
% % [jd] = date2jd(2004, 4, 6, 7, 51, 27.946047);                   % Calculate JD using JPfunction
% [r_eci,v_eci] = ecef2eci2(jd,r_ecef,v_ecef)


%% Average inertial rotation rate of the earth radians per second
omega_e = 7.292115146706979e-005;

%% find Earth orientation data
[IERS] = iers(jd);
xp = IERS(1); yp = IERS(2);
dUT1 = IERS(3); LOD = IERS(4);
ddPsi = IERS(5); ddeps = IERS(6);

%% Compute number of Julian Centuries 
[yy,mo,dd,hh,mm,ss] = jd2date2(jd);
TT = ss - 32.184;
if TT < 0
    TT = TT + 60;
    mm = mm - 1;
end
dAT = TT - 32;
if dAT < 0
    dAT = dAT + 60;
    mm = mm - 1;
end
dUT = dAT - dUT1;
if dUT < 0
    dUT = dUT + 60;
    mm = mm - 1;
end
UT1 = dAT + dUT1;
if UT1 > 60
    UT1 = UT1 - 60;
    mm = mm + 1;
end

[JD] = date2jd(yy, mo, dd, hh, mm, UT1);                        % JD_UT1 for Greenwich Apparent Sidereal Time(Theta)
Ttdb = (jd -2451545)/36525;                                     % example : p.230, vallado 4th edition
% Ttdb = 0.0426236319;                                            % example : p.230, vallado 4th edition
T2 = Ttdb^2; 
T3 = Ttdb^3; 
T4 = Ttdb^4; 

Tut1 = (JD - 2451545)/36525;                                    % Calculate Ttdb using JD_UT1's TT
Tut12 = Tut1^2; 
Tut13 = Tut1^3; 
Tut14 = Tut1^4;

%% transforamtion(Precession)
zeta  = (2306.2181/3600*Ttdb + 0.30188/3600*T2 + 0.017998/3600*T3);         % p.227, vallado 4th edition
theta = (2004.3109/3600*Ttdb - 0.42665/3600*T2 - 0.041833/3600*T3);         % p.227, vallado 4th edition
z     = (2306.2181/3600*Ttdb + 1.094/3600*T2 + 0.018203/3600*T3);           % p.227, vallado 4th edition

eci2mod = ROT3(zeta) * ROT2(-theta) * ROT3(z);                        

%% transforamtion(Nutation)
M_moon  = rem((485868.249036/3600 + 1717915923.2178/3600*Ttdb + 31.8792/3600*T2 + 0.051635/3600*T3...
    - 0.00024470/3600*T4),360);    % exampl 3-14, p. 221, reference : eq 3-82 p.225, vallado 4th edition
M_sun = rem((1287104.79305/3600 + 129596581.0481/3600*Ttdb - 0.5535/3600*T2 + 0.000136/3600*T3...
    - 0.00001149/3600*T4),360);    % exampl 3-14, p. 221, reference : eq 3-82 p.225, vallado 4th edition
uM_moon = rem((335779.526232/3600 + 1739527262.8478/3600*Ttdb - 12.7512/3600*T2 - 0.001037*T3...
    + 0.00000417*T4),360);    % exampl 3-14, p. 221, reference : eq 3-82 p.225, vallado 4th edition
D_sun = rem((1072260.70369/3600 + 1602961601.2090/3600*Ttdb - 6.3706/3600*T2 + 0.006593/3600*T3...
    - 0.00003169*T4),360);    % exampl 3-14, p. 221, reference : eq 3-82 p.225, vallado 4th edition
Omega_moon = rem((450160.398036/3600 - 6962890.5431/3600*Ttdb + 7.4722/3600*T2 + 0.007702*T3...
    - 0.00005939*T4),360);    % exampl 3-14, p. 221, reference : eq 3-82 p.225, vallado 4th edition

[Ai, Bi, Ci, Di, ai] = nutation_coef;                           % table-6, p. 1043, Vallado 4th edition
ap   = ai * [ M_moon; M_sun; uM_moon; D_sun; Omega_moon ];
dPsi = sum( (Ai*0.0001 + Bi*0.0001*Ttdb) .* sind(ap) )/3600 ;        % eq 3-83 p.226, vallado 4th edition
dEps = sum( (Ci*0.0001 + Di*0.0001*Ttdb) .* cosd(ap) )/3600;        % eq 3-83 p.226, vallado 4th edition
epsBar = (23.439291 - 0.0130042*Ttdb - 1.64e-7*T2 + 5.04e-7*T3);    % eq 3-83 p.226, vallado 4th edition

eps = epsBar + dEps ;                                                % eq 3-85 p.226, vallado 4th edition

mod2tod = ROT1(-epsBar) * ROT3(dPsi) * ROT1(eps);               % eq 3-86, p. 226, Vallado 4th edition


%% compute GAST(Greenwich Apparent Sidereal Time)
Eeq = dPsi*cosd(epsBar) +0.00264/3600*sind(Omega_moon) + 0.00063/3600*sind(2*Omega_moon);
% gmst1 = 67310.54841+(876600*60*60+8640184.812866)*Ttdb +0.093104*T2 -6.2e-6*T3;         % eq 3-47 p.188, vallado 4th edition(jd is book's(p.230) example)
gmst1 = 67310.54841+(876600*60*60+8640184.812866)*Tut1 +0.093104*Tut12 - 6.2e-6*Tut13;     % eq 3-47 p.188, vallado 4th edition(jd is jpfucntion's)
gmst2 =rem(gmst1,86400);
gmst = gmst2/240;
if gmst < 0
    gmst = gmst+360;
end

gmst12 = 24110.54841 + 8640184.812866*Tut1 +0.093104*Tut12 - 6.2e-6*Tut13;         % eq 3-45 p.188, vallado 4th edition(jd is jpfucntion's)
omegaPrec = 1.002737909350795 + 5.9006e-11*Tut1 - 5.9e-15*T2;                  % eq 3-45 p.188, vallado 4th edition(jd is jpfucntion's)
ut1 = 7*60*60+ 51*60 -(28.386009-0.4399619);
gmst_2 = rem(gmst12+omegaPrec* ut1,86400);                                      % eq 3-46 p.188, vallado 4th edition(jd is jpfucntion's)
gmst_ = gmst_2/240;
if gmst_ < 0
    gmst_ = gmst_+360;
end

gast = gmst + Eeq;                              % eq 3-79, p. 224, Vallado 4th edition
gast_ = gmst_ + Eeq;                            % eq 3-79, p. 224, Vallado 4th edition

%% transformation(Sidereal time)
tod2pef = ROT3(-gast);                          % eq 3-80, p. 224, Vallado 4th edition

%% transforamtion(Polar Motion)
xp = xp/3600; yp = yp/3600;
ecef2pef = ROT1(yp)*ROT2(xp);                   % eq 3-77, p. 223, Vallado 4th edition

om_plus= omega_e*(1-(0.0015563/86400));         % eq 3-75, p. 222, Vallado 4th edition
we=om_plus * dotR(gast);
%% Transformation
rmod = eci2mod'*r_eci';
vmod = eci2mod'*v_eci';

rtod = mod2tod'*rmod;
vtod = mod2tod'*vmod;

rpef = tod2pef'*rtod;
vpef = tod2pef'*vtod + we'*rtod;

recef = (ecef2pef' * rpef)';
vecef = (ecef2pef' * vpef)';