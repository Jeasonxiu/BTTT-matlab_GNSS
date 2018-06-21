function [estm] = PP_per_epoch_gps(Raw,data,AppPos,eleCut)
% 
%
%   function [estm] = PP_per_epoch_gps(Raw,data)
%
%   Read the data(Raw, data) and estimate Position(GPS/BDS)
%
%   input Raw : raw measurement
%   input data : eph
%   input AppPos : t-1 receiver Pos
%   input eleCut : Cut off Elevation Angle
%
%   Example : [estm] = PP_per_epoch_gps(Raw,data,AppPos,eleCut)
%
%   coded by Joonseong Gim, Sept 04, 2017
%
%

% clear all;
% close all
% %
% load('one_epoch412.mat');
% TruePos = [-3041241.741 4053944.143 3859873.640];
% AppPos = TruePos;
% eleCut = 15;
%% �Һ� ���� ����: ���� �ӵ�, ����ġ
CCC = 299792458.;   % CCC = Speed of Light [m/s]
%% �Ӱ���� ����
%% ���� �ڵ鸵
LeapSecBDS = 14;

%% ����� ����ġ ����
QMfile = Raw{3,1}(find(Raw{3,1}(:,1) ~= 0),:);
for i=1:length(QMfile(:,1))
    prn = QMfile(i,2);
    if prn < 200
        TYPE = 120;
    elseif prn > 200
        TYPE = 220;
    end
    QM(i,:)=[QMfile(i,1), prn, TYPE, QMfile(i,3)];
    QMsnr(i,:)=[QMfile(i,1), prn, TYPE, QMfile(i,7)];
end
%% Eph ����
eph = data{3,1}(find(data{3,1}(:,28) == 1 & data{3,1}(:,29) == 1 & data{3,1}(:,30) == 1 & data{3,1}(:,18) ~= 204),:);
%% Klobucher �� ����
al = data{2}(1:4);
be = data{2}(5:8);

%% GPS week, GPS second ����
gw = Raw{2}(2);
gs = Raw{2}(1);
%% �ʱ���ǥ ȹ��
gd = xyz2gd(AppPos); AppLat = gd(1); AppLon = gd(2);
x = AppPos;

%% SV pos
PRN1e = intersect(QM(:,2), eph(:,18));
NoSats = length(PRN1e);
for i = 1:NoSats
    prn = PRN1e(i);
    icol = PickEPH(eph, prn, gs);
    STT = GetSTTbrdm(gs, eph, icol, x(1:3)'); % ��ȣ���޽ð� ���
    
    if prn > 100 && prn < 200
        tc = gs - STT ;             % ��ȣ���޽ð� ����
    elseif prn > 200 && prn < 300
        tc = gs - STT- LeapSecBDS;  % ��ȣ���� �ð� ����... bds���� �ݿ�
    end
    
    SatPos = GetSatPosNC_GC(eph, icol, tc); % ������ġ ����
    SatPos = RotSatPos(SatPos, STT); % ��������ȿ�� ���
    
    vec_rho = SatPos - x(1:3)';
    rho = norm(vec_rho);
    [az,el] = xyz2azel(vec_rho, AppLat, AppLon);
    SV(i,:) = [prn, SatPos', az, el];
end
% SV = SV(find(SV(:,1) < 200),:);
% SV = SV(find(SV(:,1) ~= 209),:);

%% ���� ������ Ȯ���Ͽ� GPS/BDS or GPS or BDS ���� Ȯ��
NumGPSSV = length(SV(find(SV(:,1) < 200 & SV(:,6) > eleCut)));
NumBDSSV = length(SV(find(SV(:,1) > 200 & SV(:,6) > eleCut)));
SV = SV(find(SV(:,6) > eleCut),:);
%% ������ �ʿ��� �ʱ�ġ ����
Maxiter = 10;
EpsStop = 1e-5;
ctr = 1; ctr2 = 1;
NumGPSsUsed = 0;
NumBDSsUsed = 0;
NumGLOsUsed = 0;
%% �������� ����
estm = [];
tic;

if NumGPSSV >= 1 & NumBDSSV == 0
    x = [x, ctr];    H = zeros(1,4);
    GPSflag = 1; BDSflag = 0;
elseif NumGPSSV == 0 & NumBDSSV >= 1
    x = [x, ctr];    H = zeros(1,4);
    GPSflag = 0; BDSflag = 1;
elseif NumGPSSV >= 1 & NumBDSSV >= 1
    x = [x, ctr, ctr2];    H = zeros(1,5);
    GPSflag = 1; BDSflag = 1;
end

for iter = 1:Maxiter
    %% ���� ���� 4�� �����϶� ��� ����
    if length(SV(:,1)) <= 4
        break;
    end
    vec_site = x(1:3)';
    ZHD = TropGPTh(AppPos, gw, gs); %: TROP: GPT
    NumGPSsUsed = 0;
    NumBDSsUsed = 0;
    NumGLOsUsed = 0;
    
    cnt2=1;
    
    for i=1:length(SV(:,1))
        prn = SV(i,1);
        obs = QM(find(QM(:,2) == prn),4);
        icol = PickEPH(eph, prn, gs);
        toe = eph(icol, 8); a = eph(icol, 19); b = eph(icol, 20); c = eph(icol, 21); Tgd = eph(icol, 23); Toc = eph(icol, 26);
%         SatPos = SV(find(SV(:,1) == prn),2:4);
        STT = GetSTTbrdm(gs, eph, icol, x(1:3)'); % ��ȣ���޽ð� ���
        
        if prn > 100 && prn < 200
            tc = gs - STT ;             % ��ȣ���޽ð� ����
        elseif prn > 200 && prn < 300
            tc = gs - STT- LeapSecBDS;  % ��ȣ���� �ð� ����... bds���� �ݿ�
        end
        SatPos = GetSatPosNC_GC(eph, icol, tc); % ������ġ ����
        SatPos = RotSatPos(SatPos, STT); % ��������ȿ�� ���
        
        vec_rho = SatPos' - vec_site';
        [az,el] = xyz2azel(vec_rho, AppLat, AppLon);
        rho = norm(vec_rho);
        if prn < 200 | (prn > 200 & prn < 300)
            dIono = klobuchar(al, be, gs, az, el, vec_site); % �̿��� ����(Klobuchar ��)
            dTrop = ZHD2SHD(gw, gs, AppPos, el, ZHD); % ����� ����
        else
        end
        dRel = GetRelBRDC(eph, icol, tc); % ��뼺ȿ��
        dtSat = a + b*(tc - Toc) + c*(tc - Toc)^2 - Tgd + dRel; % �����ð���� ���... �׷������, ��뼺ȿ�� ����
        W(cnt2,cnt2) = 1;
        if (GPSflag == 1 & BDSflag == 0) | (GPSflag == 0 & BDSflag == 1)
            com = rho + x(4) - CCC*dtSat + dTrop + dIono; % GPS or BDS ��갪
            H(cnt2,:) = [-vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 1];
            if prn < 200
                NumGPSsUsed = NumGPSsUsed + 1;
            elseif prn > 200 & prn < 300
                NumBDSsUsed = NumBDSsUsed + 1;
            end
            Constellation = 1;
        elseif GPSflag == 1 & BDSflag ==1
            if prn < 200
                com = rho + x(4) - CCC*dtSat + dTrop + dIono; % GPS ��갪
                H(cnt2,:) = [-vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 1 0];
                NumGPSsUsed = NumGPSsUsed + 1;
            elseif prn > 200 & prn < 300
                com = rho + x(5) - CCC*dtSat + dTrop + dIono; % GPS ��갪
                H(cnt2,:) = [-vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 0 1];
                NumBDSsUsed = NumBDSsUsed + 1;
            else
                
            end
            Constellation = 2;
        end
        y(cnt2,1) = obs - com;
        cnt2=cnt2+1;
        
    end
    HTH = H(1:cnt2-1,:)'*inv(W(1:cnt2-1,1:cnt2-1))*H(1:cnt2-1,:);
    HTy = H(1:cnt2-1,:)'*inv(W(1:cnt2-1,1:cnt2-1))*y(1:cnt2-1,:);
    xhat = inv(HTH) * HTy;
    x = x + xhat';
    if norm(xhat) < EpsStop;
        estm(1,1) = gs;
        if Constellation < 2
            estm(1,2:5) = x(1:4);
        else
            estm(1,2:6) = x(1:5);
        end
        estm(1,7) = NumGPSsUsed;
        estm(1,8) = NumBDSsUsed;
        estm(1,9) = NumGLOsUsed;
        fprintf('%8d : %3d : %8.2f : %8.2f : %8.2f\n', j, iter, x(1)' - AppPos(1), x(2)' - AppPos(2),x(3)' - AppPos(3));
        break;
    end
end


toc;
if isempty(estm)
    estm(1,1) = gs;
    estm(1,2:4) = AppPos;
    estm(1,7) = NumGPSsUsed;
    estm(1,8) = NumBDSsUsed;
    estm(1,9) = NumGLOsUsed;
end
% [dXYZ, dNEV] = PosTErrors2(estm(:,1), TruePos, estm(:,2:5));
