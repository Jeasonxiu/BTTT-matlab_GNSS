function [Trop] = Hopfield(el, C, vec_site, p)
% function [Trop] = Hopfield(el, C, vec_site, p);
%
%   Troposphere correction using Hopfield Model
%   
%   input el : Elevation angle of GPS Satellite(degree)
%   input T : Temperature(Celsius)
%   input vec_site : Receiver XYZ(ECEF)
%   input p : Pressure(if not exist p = 9999)
%
%   Example : [Trpo] = Hopfield(el, 23, vec_site, 9999);
%
%   coded by Joonseong Gim, Feb 11, 2016
%
e = 0;
[h] = xyz2gd(vec_site);             % h계산
T = C + 273.16;                     % Celsius to Kelvin
D2R = pi / 180;                     % degree to radian

if p == 9999
    p = (-0.12*h(3))+1027.6;        % 입력된 대기압이 없을 경우 한반도형 대기모델 이용
end

hd = 40136 +148.72-(T-273.16);
hw = 11000;
mapFd = (sin(sqrt(el^2+6.25)*D2R))^(-1);
mapFw = (sin(sqrt(el^2+2.25)*D2R))^(-1);
Trop_d = (10^(-6)/5)*(77.64*p/T)*mapFd*hd;
Trop_w = (10^(-6)/5)*(-12.96*T+3.718*10^5)*mapFw*(e/T^2)*hw;
Trop = Trop_d + Trop_w;