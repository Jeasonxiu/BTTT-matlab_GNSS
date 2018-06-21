function [Trop] = MHopfield(el, C, vec_site, p)
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
[h] = xyz2gd(vec_site);             % h계산
e = 0;
if p == 9999
    p = (-0.12*h(3))+1027.6;        % 입력된 대기압이 없을 경우 한반도형 대기모델 이용
end

T = C + 273.16;                     % Celsius to Kelvin
RE = 6378.137                       % Goad & Goodman 1974 참조
hd = 40136+148.72*(T-273.16);
hw = 11000;

%% hydrostatic, wet parameters
rid = sqrt((RE+hd)^2-(RE*cosd(el)^2));
riw = sqrt((RE+hw)^2-(RE*cosd(el)^2));
ad = -(sind(el)/hd); bd = -((cosd(el)^2)/(2*hd*RE)));
aw = -(sind(el)/hw); bw = -((cosd(el)^2)/(2*hw*RE)));

alphad = [1 4*ad 6*ad^2+4*bd 4*ad*(ad^2+3*bd) ad^4+12*ad^2*bd+6*bd^2...
    4*ad*bd*(ad^2+3*bd) bd^2*(6*ad^2+4*bd) 4*ad*bd^3 bd^4];
alphaw = [1 4*aw 6*aw^2+4*bw 4*aw*(aw^2+3*bw) aw^4+12*aw^2*bw+6*bw^2...
    4*aw*bw*(aw^2+3*bw) bw^2*(6*aw^2+4*bw) 4*aw*bw^3 bw^4];

ALD = 0; ALW = 0;
for i = 1:9
    ALD = ALD + alphad(i);
    ALW = ALW + alphaw(i);
end


Trop_d = (10^(-6)/5)*((77.64*p/T)/sind(sqrt(el^2)+6.25))*(40136+148.72*(T-273.16));
Trop_w = (10^(-6)/5)*((-12.96*T+3.718*10^5)/sind(sqrt(el^2+2.25)))*(e/T^2)*11000;
Trop = Trop_d + Trop_w;