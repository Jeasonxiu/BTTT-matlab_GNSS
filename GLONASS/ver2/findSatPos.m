function [epochSat,PoS,Vel] = findSatPos(Sat_ar,tc,prn)

% function [PoS,Vel] = findSatPos(Sat_ar,tc,prn)
%
% input : GatSatPosGLO_new�� 1�� �������� ������� Sat array
%         prn : Satellite number
%         tc : gs - LeapSecend
%
% April 6, 2015, Mi-So Kim

%   Sat_ar = [1prn 2�ð� 3x 4y 5z 6vx 7vy 8vz 9LSx 10LSy 11LSz]

% �ش� ������ȣ ã��
SatP1 = find(Sat_ar(:,1) == prn);

% �ش� �ð� ã��
% �� �ð��� ���ؼ� ã�� �ð��� ���� �� �� �ð��� ��, ���� �� ���� ����� �ð� ��

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