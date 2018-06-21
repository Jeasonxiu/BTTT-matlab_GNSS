function [Trop] = Corr_Trop_Hopfield(P_rec, P_sat, P)

% Hopfield Model(1969)
% 지표면의 굴절지수를 이용하여 고도에 따른 굴절률 변화를 나타낸 모델
% 입력변수: 관측소 좌표 P_rec(m), 위성 좌표 P_sat(m). 대기압 P(hPa)
%                대기압 값이 없을 경우 9999입력
% 출력변수: 대류권 신호지연량

T = 273.16;
e = 0;

% 관측소의 3차원 좌표 P_rec를 위도, 경도, 타원체고로 변환
[lat, lon, h] = xyz2gd(P_rec);

% 대기압 입력값이 없으면 한반도형 대기모델 사용
if P == 9999
    P = (-0.12*h)+1027.6;        % 입력된 대기압이 없을 경우 한반도형 대기모델 이용
end

% 고도각 E 산출
delta = P_sat - P_rec;
[topo] = xyz2topo(delta, lat, lon);
[AzEl] = topo2AzEl(topo);
E = AzEl(2);

Trop_d = (10^(-6)/5)*(77.64/(sind(sqrt(E^2+6.25))))*P/T*(40136+148.72*(T-273.16));
Trop_w = (10^(-6)/5)*((-12.96*T + 3.718*10^5)/(sind(sqrt(E^2+2.25)))) * (e/(T^2)) * 11000;
Trop = Trop_d + Trop_w;