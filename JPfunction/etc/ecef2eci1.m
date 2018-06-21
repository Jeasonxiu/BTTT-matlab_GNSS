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
%   coded by Joonseong Gim, April 4, 2016
%



% Compute Sidereal Angle
% theta = jd2STA(jd);
theta = JD2GAST(jd);
% Average inertial rotation rate of the earth radians per second
omega_e = 7.29211585275553e-005;
% Compute r_eci Components
cost=cosd(theta);sint=sind(theta);
r_eci(:,1)=cost.*r_ecef(:,1)-sint.*r_ecef(:,2);
r_eci(:,2)=sint.*r_ecef(:,1)+cost.*r_ecef(:,2);
r_eci(:,3)=r_ecef(:,3);
% Compute v_eci Components
v_eci(:,1)=cost.*v_ecef(:,1)-sint.*v_ecef(:,2)-omega_e*sint.*r_ecef(:,1)-omega_e*cost.*r_ecef(:,2);
v_eci(:,2)=sint.*v_ecef(:,1)+cost.*v_ecef(:,2)+omega_e*cost.*r_ecef(:,1)-omega_e*sint.*r_ecef(:,2);
v_eci(:,3)=v_ecef(:,3);