function [Trop] = Corr_Trop_Saastamoinen(P_rec, P_sat, P)

% Saastamoinen Model(1973)
% 고도각에 따른 대류층 굴절률을 추정한 모델
% 입력변수: 관측소 좌표 P_rec(m), 위성 좌표 P_sat(m). 대기압 P(hPa)
%                대기압 값이 없을 경우 9999입력
% 출력변수: 대류층 신호지연량(m)

T = 273.16;
e = 0;

% 관측소의 3차원 좌표 P_rec를 위도, 경도, 타원체고로 변환
[lat, lon, h] = xyz2gd(P_rec);

% 대기압 입력값이 없으면 한반도형 대기모델 사용
if P== 9999
    P = (-0.12*h)+1027.6;         % 입력된 대기압 값이 없을 경우 한반도형 대기모델 이용
end

% 고도각 E 산출
delta = P_sat - P_rec;
[topo] = xyz2topo(delta, lat, lon);
[AzEl] = topo2AzEl(topo);
E = AzEl(2);

Z = 90-E;

Trop = (0.002277/cosd(Z))*(P+(1255/T+0.05)*e-tand(Z)^2);