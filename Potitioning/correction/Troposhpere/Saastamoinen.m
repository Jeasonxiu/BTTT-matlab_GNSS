function [Trop] = Saastamoinen(el, C, vec_site, p)
% function [Trop] = Saastamoinen(el, C, vec_site, p)
%
%   Troposphere correction using Saastamoinen Model
%   
%   input el : Elevation angle of GPS Satellite(degree)
%   input T : Temperature(Celsius)
%   input vec_site : Receiver XYZ(ECEF)
%   input p : Pressure(if not exist p = 9999)
%
%   Example : [Trop] = Saastamoinen(el, 2, vec_site, 9999);
%
%   coded by Joonseong Gim, Feb 11, 2016
%
e = 0;
T = C + 273.16;                     % Celsius to Kelvin
[h] = xyz2gd(vec_site);             % h계산
D2R = pi/180;
if p == 9999
    p = (-0.12*h(3))+1027.6;        % 입력된 대기압이 없을 경우 한반도형 대기모델 이용
end

Z = (90 - el) * D2R;

Trop = (0.002277/cos(Z))*(p+(1255/T+0.05)*e-tan(Z)^2);
