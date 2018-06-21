function [epochSat,PoS,Vel] = findSatPos(Sat_ar,tc,prn)

% function [PoS,Vel] = findSatPos(Sat_ar,tc,prn)
%
% input : GatSatPosGLO_new의 1초 간격으로 만들어진 Sat array
%         prn : Satellite number
%         tc : gs - LeapSecend
%
% April 6, 2015, Mi-So Kim

%   Sat_ar = [1prn 2시간 3x 4y 5z 6vx 7vy 8vz 9LSx 10LSy 11LSz]

% 해당 위성번호 찾기
SatP1 = find(Sat_ar(:,1) == prn);

% 해당 시간 찾기
% → 시간을 비교해서 찾는 시간이 있을 땐 그 시간을 겟, 없을 땐 가장 가까운 시간 겟

n = length(SatP1);
if n == 0
    return;
end

tk = SatP1(1);
dtmin = Sat_ar(tk,2) - tc;
for k1 = 2 : n
    kk = SatP1(k1);
    dt = Sat_ar(kk,2) - tc;
    %     if dt < 0
    if dt <= 0
        tk = kk;
        dtmin = dt;
    end
    % end
end

epochSat = Sat_ar(tk,:);
PoS = epochSat(3:5); Vel = epochSat(6:8);