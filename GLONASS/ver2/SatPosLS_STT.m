function [SatPos, SatVel] = SatPosLS_STT(Sat_ar,gs,prn,LeapSecond,STT,TauC)

% <input>
%    epochSat : function-findSatPos����
%    Leapsecond : ����
%    STT : ���� ��ȣ���� �ð�
%
% <output>
%    SatPos, SatVel : ���ʿ� ��ȣ���� �ð��� ������ ������ ��ġ�� �ӵ�
%
%  Code by Miso Kim, April 8th, 2015

tc = gs - LeapSecond;

[epoch,SatPos,SatVel] = findSatPos(Sat_ar,tc,prn);


if epoch(2) == tc % Sat_ar�� gs - LeapSecond�� �ִ� ���: STT�� ����
    deltat = - STT + TauC;
%     [xyz,vxyz,xyzLS] = RK4_new(epoch,deltat);
    [xyz,vxyz,xyzLS] = RK4_my(epoch,deltat);
    SatPos = xyz; SatVel = vxyz;
else % Sat_ar�� gs - LeapSecond�� ���� ���: gs�� tc�� ���̸�ŭ(LS), STT ����
    k = 0;
    intv = epoch(2) - tc;
    x = epoch;
    if intv > 0
        deltat = -1;
    else
        deltat = 1;
    end
    for i = 1:intv
        k = k + 1;
%         [xyz,vxyz,xyzLS] = RK4_new(x,deltat);
        [xyz,vxyz,xyzLS] = RK4_my(x,deltat);
        SatR(k,1) = x(1); SatR(k,2) = x(2) + deltat;
        SatR(k,3:5) = xyz; SatR(k,6:8) = vxyz; SatR(k,9:11) = xyzLS;
        x = SatR(k,:);
    end
    deltat = - STT + TauC;
%     [xyz,vxyz,xyzLS] = RK4_new(SatR(k,:),deltat);
    [xyz,vxyz,xyzLS] = RK4_my(SatR(k,:),deltat);
    SatPos = xyz; SatVel = vxyz;
end