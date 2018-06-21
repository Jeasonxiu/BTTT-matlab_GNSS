function [estm] = PP_per_epoch(Raw,data,AppPos,eleCut)
%
%
%   function [estm] = PP_per_epoch(Raw,data)
%
%   Read the data(Raw, data) and estimate Position(GPS/BDS)
%
%   input Raw : raw measurement
%   input data : eph
%
%   Example : [estm] = PP_per_epoch(Raw,data,AppPos)
%
%   coded by Joonseong Gim, Sept 04, 2017
%
%

% clear all;
% close all
%
% load('one_epoch.mat');
% TruePos = [-3041241.741 4053944.143 3859873.640];
% AppPos = TruePos;
%% �Ӱ���� ����
% eleCut = 0;
%% �Һ� ���� ����: ���� �ӵ�, ����ġ
CCC = 299792458.;   % CCC = Speed of Light [m/s]
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

%% GPS week ����
gw = Raw{2}(2);

%% �ʱ���ǥ ȹ��
gd = xyz2gd(AppPos); AppLat = gd(1); AppLon = gd(2);
x = AppPos;

%% ������ �ʿ��� �ʱ�ġ ����
Maxiter = 10;
EpsStop = 1e-5;
ctr = 1; ctr2 = 1;


%% �������� ����
nEst = 0;
estm = [];
tic;
%% Leat Square
PRN1e = intersect(QM(:,2), eph(:,18));
NoSats = length(PRN1e);
NoGPSSV = length(find(PRN1e(:) < 200));
NoBDSSV = length(find(PRN1e(:) > 200));
if NoGPSSV ~= 0 & NoBDSSV ~= 0
    H = zeros(1,5);
    x = [AppPos 1 1]';
else
    H = zeros(1,4);
    x = [AppPos 1]';
end
for iter = 1:Maxiter
    %% Raw measurement �� Ȯ���Ͽ� GPS/BDS or GPS or BDS ���� Ȯ��
    
    gs = QM(1,1);
    if NoSats <= 4
        break
    end
    vec_site = x(1:3)';
    ZHD = TropGPTh(AppPos, gw, gs); %: TROP: GPT
    NumGPSsUsed = 0;
    NumBDSsUsed = 0;
    NumGLOsUsed = 0;
    
    cnt2=1;
    for i = 1:NoSats
        
        prn = PRN1e(i);
        obs = QM(find(QM(:,2) == prn),4);
        snr = QMsnr(find(QMsnr(:,2) == prn),4);
        icol = PickEPH(eph, prn, gs);
        toe = eph(icol, 8); a = eph(icol, 19); b = eph(icol, 20); c = eph(icol, 21); Tgd = eph(icol, 23); Toc = eph(icol, 26);
        STT = GetSTTbrdm(gs, eph, icol, x(1:3)); % ��ȣ���޽ð� ���
        
        if prn > 100 && prn < 200
            tc = gs - STT ;             % ��ȣ���޽ð� ����
        elseif prn > 200 && prn < 300
            tc = gs - STT- LeapSecBDS;  % ��ȣ���� �ð� ����... bds���� �ݿ�
        end
        dRel = GetRelBRDC(eph, icol, tc); % ��뼺ȿ��
        dtSat = a + b*(tc - Toc) + c*(tc - Toc)^2 - Tgd + dRel; % �����ð���� ���... �׷������, ��뼺ȿ�� ����
        SatPos = GetSatPosNC_GC(eph, icol, tc); % ������ġ ����
        SatPos = RotSatPos(SatPos, STT); % ��������ȿ�� ���
        
        vec_rho = SatPos - vec_site';
        rho = norm(vec_rho);
        [az,el] = xyz2azel(vec_rho, AppLat, AppLon);
        
        if prn > 100 && prn < 200
            dIono = klobuchar(al, be, gs, az, el, x(1:3)); % �̿��� ����(Klobuchar ��)
            dTrop = ZHD2SHD(gw, gs, AppPos, el, ZHD); % ����� ����
        elseif prn > 200 && prn < 300
            dIono = klobuchar(al, be, gs, az, el, x(1:3)); % �̿��� ����(Klobuchar ��)
            dTrop = ZHD2SHD(gw, gs, AppPos, el, ZHD); % ����� ����
        end
        %             dIono = 0;
        %             dTrop = 0;
        if el >=eleCut % �Ӱ����
            W(cnt2,cnt2) = 1;
            %                 W = MakeW_elsnr(el,snr);
            if prn > 100 && prn < 200
                com = rho + x(4) - CCC*dtSat + dTrop + dIono; % gps ��갪
                if length(H(1,:)) == 5
                    H(cnt2,:) = [-vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 1 0];
                else
                    H(cnt2,:) = [-vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 1];
                end
                NumGPSsUsed = NumGPSsUsed + 1;
            elseif prn > 200 && prn < 300
                if length(H(1,:)) == 5
                    com = rho + x(5) - CCC*dtSat + dTrop + dIono; % bds ��갪
                    H(cnt2,:) = [-vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 0 1];
                else
                    com = rho + x(4) - CCC*dtSat + dTrop + dIono; % bds ��갪
                    H(cnt2,:) = [-vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 1];
                end
                NumBDSsUsed = NumBDSsUsed + 1;
            end
            y(cnt2,1) = obs - com;
            cnt2=cnt2+1;
            
        end
    end
    HTH = H(1:cnt2-1,:)'*inv(W(1:cnt2-1,1:cnt2-1))*H(1:cnt2-1,:);
    HTy = H(1:cnt2-1,:)'*inv(W(1:cnt2-1,1:cnt2-1))*y(1:cnt2-1,:);
    xhat = inv(HTH) * HTy;
    x = x + xhat;
    if norm(xhat) < EpsStop;
        %     if iter == Maxiter
        nEst = nEst + 1;
        estm(nEst,1) = gs;
        estm(nEst,2:5) = x(1:4);
        estm(nEst,7) = NumGPSsUsed;
        estm(nEst,8) = NumBDSsUsed;
        estm(nEst,9) = NumGLOsUsed;
        %             Scount(nEst,1) = NoUsed_gps;
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
