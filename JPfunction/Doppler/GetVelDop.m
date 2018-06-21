function [VelDop] = GetVelDop(gs, prns, dops, eph, vec_site)
%
%function [VelDop] = GetVelDop(gs, prns, dops, eph, vec_site)
%
% DO: Compute velocity based on Doppler measurements
%
% <input>   gs: GPS week second
%           prns: array of prns observed at the epoch
%           dops: doppler measurements at the epoch correspong to the prn
%           eph: ephemeris array read using ReadEPH
%           vec_site: receiver position at the epoch
%
% <output>  VelDop: Velocity computed from Doppler observables
%
% Copyright: Kwan-Dong Park@Jipyong Space, 7/12/2015
%
%% 필요한 상수 및 매개변수 정의
CCC = 299792458.; % CCC = Speed of Light [m/s]
f1 = 1575.42e6;   % L1 주파수
lmbd1 = CCC/f1;   % L1 파장
%% 초기값 4개를 모두 0으로 설정해도 큰 문제 없는 것으로 파악하였음
x = [0 0 0 0]; x = x';
vel_site = x(1:3)';
%% 정규방정식에 필요한 행렬 초기화
HTH = zeros(4,4);
HTy = zeros(4,1);
%% 선형 최소제곱법 추정 과정
NoSats = length(prns);
for kS = 1:NoSats
    prn = prns(kS);
    obs = -dops(kS);
    
    icol = PickEPH(eph, prn, gs);
    [vec_sat,vel_sat] = GetSatVelNC(eph, icol, gs);
    
    vec_site2sat = vec_sat - vec_site;
    vel_site2sat = vel_sat - vel_site;
    
    com = dot(vec_site2sat/norm(vec_site2sat), vel_site2sat) + x(4);
    y = lmbd1*obs - com;
    
    H = [ -vec_site2sat/norm(vec_site2sat) 1];
    
    HTH = HTH + H'*H;
    HTy = HTy + H'*y;
end
%% 선형문제이므로 반복 계산 필요 없음
x = x + inv(HTH)*HTy;
VelDop = x(1:3)';
