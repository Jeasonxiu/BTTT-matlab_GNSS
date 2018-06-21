function [theta] = jd2STA(JD)
% 
%   input jd[1 x 1] : Julian day
%
%   output theta : Sidereal Time Angle(degree)
%
%   coded by Joonseong Gim, April 6, 2016
%   

% Compute T
% t=(jd-2451545)/36525;
% % Compute Sidereal Angle
% theta=(280.46061837+360.98564736629*(jd-2451545)+0.000387933*t^2-t^3/38710000);

[yy,mo,dd,hh,mm,ss] = jd2date(JD);
timehack = [yy,mo,dd,hh,mm,ss];
D = JD - 2451545;

% Calculate UT in decimal hours
UT=timehack(1,4)+(timehack(1,5)+timehack(1,6)/60)/60;
% Calculate T (number of centuries since UT began)
T=(JD-2451545)/36525;
% Calculate T0
T0 = 6.697374558 + (2400.051336*T) + (0.000025862*T^2) + (UT*1.0027379093);

% Reduce T0 to a value between 0 and 24 by adding or subtracting multiples of 24
T0=mod(T0,24);
% Greenwich mean sidereal time in decimal hours
GMST=T0;

% Solve for Apparent Sidereal time by solving for the
% nutation in right ascension
Om=125.04-0.052954*D;                              %Longitude of the ascending node of the Moon
L=280.47+0.98565*D;                                %Mean Longitude of the Sun
eps=23.4393-0.0000004*D;                           %obliquity
delta_psi=-0.000319*sind(Om)-0.000024*sind(2*L);   %nutation in longitude
eqeq=delta_psi*cosd(eps);                          %equation of the equinoxes
GAST=GMST + eqeq;                                  %Greenwich Apparent Sidereal Time

% Convert from decimal hours to degs to rads to match input for conversion
% to/from Earth Centered Inertial/Earth Centered Rotational functions
GAST=GAST*360/24;
theta = GAST;