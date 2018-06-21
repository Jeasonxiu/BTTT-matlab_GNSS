function [KE] = rv2KE(JD,r,v)
%
% function [KE] = rv2KE(JD,R,V)
%
%   Determining Keplerian Elements from GLONASS Satellites's POS, VEL
%   
%   input JD[1 x 1] : Julian day for ECEF2ECI function
%   input R[1 x 3] : Position Vector in ECEF(PZ-90) coordinate frame of reference
%   input V[1 x 3] : Velocity vector in ECEF(PZ-90) coordinate frame of reference
%
%   output KE[a,e,i,w,O,f,p] Kepliarn Elements Matrix
%   i,w,O,f(radians)
%
%   Example : [KE] = rv2KE(JD,R,V)
%
%   coded by Joonseong Gim, April 4, 2016
%

%% 함수 돌리기 전 테스트를 위한 파라미터 생성
% clear all
% close all
% load('glomat.mat');
% R = EphGlo(1,3:5);                                              % Position Vector
% V = EphGlo(1,6:8);                                              % Velocity Vector
% JD = EphGlo(1,18);                                              % Velocity Vector
% r = [6524834, 6862875, 6448296];
% v = [4901.327, 5533.756, -1976.341];
%% r,v(PZ-90) to R,V(WGS84)
R = PZ2WGS(r);                                             % Position Vector(WGS84)
V = PZ2WGS(v);                                             % Velocity Vector(WGS84)
%% r,v(WGS84) to ECI
[r v] = ecef2eci2(JD,R,V);                                       % convert ecef to eci(Joon's code)
%% Determining Keplian Elements 
mu = 3.986004418e14;                                              % Earth's gravitational constant
h = cross(r,v);                                                 % obital momentum
k=[0,0,1];
n = cross(k,h);                                                 %  
nMag = norm(n);                                                 % magnitude of the n vector                        
vMag = norm(v);                                                 % magnitude of the velocity vector                        
rMag = norm(r);                                                 % magnitude of the position vector                        
hMag = norm(h);                                                 % magnitude of the obital momentum vector 

evec = (1/mu)*((vMag^2 - mu/rMag) * r - dot(r,v) * v);          % eccentricity vector
eMag = norm(evec);                                              % magnitude of the eccentricity vector                        

EE = (vMag^2)/2 - mu/rMag;       % Energy Equation
% a = -mu/(2*EE);                                                 % Semi-Major Axis(m)
a = hMag^2/((1-eMag^2)*mu);                                     % Semi-Major Axis(m)

%Compute the angles
i = acos(h(3)/hMag);                                            % inclination(radians)
O = acos(n(1)/nMag);                                            % ascending node(radians)
if n(2) < 0
    O = 2*pi - O;
end
w = acos(dot(n,evec)/(nMag*eMag));                              % Argument of perigee(radians)
if evec(3) < 0
    w = 2*pi - w;
end
f = acos(dot(evec,r)/(eMag*rMag));                              % true Anormaly(radians)
if  dot(r,v) < 0
    f = 2*pi - f;
end
p = hMag^2/mu;
KE = [a, eMag, i, w, O, f, p];