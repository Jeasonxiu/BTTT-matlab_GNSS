%% GPS & BDS ���� ���� �ڵ� (DGNSS)
tic
clear all; close all;
warning off;
%% ��� ����
CCC = 299792458.;   % CCC = Speed of Light [m/s]
%% QM ����, ��¥ �ڵ鸵
%-----------161220 �׽�Ʈ---------%
% RTCM=load('161220_PPS1_PRC');
% DOY = 355;YY  = 16;
% FileQM='A_Pole_a';
% FileQM='A_Pole_b';
% 
% FileQM='A_��_a';
% FileQM='A_��_b';
% 
% FileQM='A_��_a';
% FileQM='A_��_b';

% FileQM='A_��_a';
% FileQM='A_��_b';

% FileQM='A_��_a';
% FileQM='A_��_b';

% 
% FileQM='B_Pole_a';
% FileQM='B_Pole_b';
% 
% FileQM='B_��_a';
% FileQM='B_��_b';

% FileQM='B_��_a';
% FileQM='B_��_b';

% FileQM='B_��_a';
% FileQM='B_��_b';

% FileQM='B_��_a';
% FileQM='B_��_b';

%-----------161223 �׽�Ʈ---------%
% load('161223_PRC.mat');
% DOY = 358;YY  = 16;
% FileQM='test_p1_L';
% FileQM='test_p1_R';

% FileQM='test_p2_L';
% FileQM='test_p2_R';
% 
% FileQM='test_p3_L';
% FileQM='test_p3_R';
% 
% FileQM='test_p4_L';
% FileQM='test_p4_R';
% 
% FileQM='test_j1_L';
% FileQM='test_j1_R';
% 
% FileQM='test_j2_L';
% FileQM='test_j2_R';
% 
% FileQM='test_j3_L';
% FileQM='test_j3_R';
% 
% FileQM='test_j4_L';
% FileQM='test_j4_R';

%-----------161227 �׽�Ʈ---------%
% RTCM = load('PPS1_161227.t41');
% DOY = 362;YY  = 16;
% FileQM='A_P_L';
% FileQM='A_P_R';

% FileQM='A_N_L';
% FileQM='A_N_R';
% 
% FileQM='A_E_L';
% FileQM='A_E_R';
% 
% FileQM='A_S_L';
% FileQM='A_S_R';
% 
% FileQM='A_W_L';
% FileQM='A_W_R';
% 
% FileQM='B_P_L';
% FileQM='B_P_R';

% FileQM='B_N_L';
% FileQM='B_N_R';
% 
% FileQM='B_E_L';
% FileQM='B_E_R';
% 
% FileQM='B_S_L';
% FileQM='B_S_R';
% 
% FileQM='B_W_L';
% FileQM='B_W_R';

%-----------161228 �׽�Ʈ---------%
% RTCM = load('PPS1_161228.t41');
% DOY = 363;YY  = 16;
% FileQM='C_P_L';
% FileQM='C_P_R';

% FileQM='C_N_L';
% FileQM='C_N_R';
% 
% FileQM='C_E_L';
% FileQM='C_E_R';
% 
% FileQM='C_S_L';
% FileQM='C_S_R';
% 
% FileQM='C_W_L';
% FileQM='C_W_R';
% 
% FileQM='D_P_L';
% FileQM='D_P_R';

% FileQM='D_N_L';
% FileQM='D_N_R';
% 
% FileQM='D_E_L';
% FileQM='D_E_R';
% 
% FileQM='D_S_L';
% FileQM='D_S_R';
% 
% FileQM='D_W_L';
% FileQM='D_W_R';

%-----------170112 �׽�Ʈ---------%
% RTCM = load('PPS1_170112.t41');
% DOY = 012;YY  = 17;
% FileQM='P_10Hz_L1';
% FileQM='P_10Hz_R1';

% FileQM='P_10Hz_L2';
% FileQM='P_10Hz_R2';

% FileQM='P_10Hz_L3';
% FileQM='P_10Hz_R3';

% FileQM='P_10Hz_L4';
% FileQM='P_10Hz_R4';

% FileQM='J_1Hz_L1';
% FileQM='J_1Hz_R1';

% FileQM='J_1Hz_L2';
% FileQM='J_1Hz_R2';

% FileQM='J_1Hz_L3';
% FileQM='J_1Hz_R3';

% FileQM='J_1Hz_L4';
% FileQM='J_1Hz_R4';

% FileQM='J_5Hz_L1';
% FileQM='J_5Hz_R1';

% FileQM='J_5Hz_L2';
% FileQM='J_5Hz_R2';

% FileQM='J_5Hz_L3';
% FileQM='J_5Hz_R3';

% FileQM='J_5Hz_L4';
% FileQM='J_5Hz_R4';

% FileQM='J_10Hz_L1';
% FileQM='J_10Hz_R1';

% FileQM='J_10Hz_L2';
% FileQM='J_10Hz_R2';

% FileQM='J_10Hz_L3';
% FileQM='J_10Hz_R3';

% FileQM='J_10Hz_L4';
% FileQM='J_10Hz_R4';
%-----------170116 �׽�Ʈ---------%
RTCM = load('PPS1_170116.t41');
DOY = 016;YY  = 17;
%----���� Type A----%
% FileQM='test1_L';
% FileQM='test1_R';

%----���� Type B----%
% FileQM='test2-1_L';
% FileQM='test2-1_R';

% FileQM='test2-2_L';
% FileQM='test2-2_R';

%----SUV Type A----%
FileQM='test3_L';
% FileQM='test3_R';

%----SUV Type B----%
% FileQM='test4-1_L';
% FileQM='test4-1_R';

% FileQM='test4-2_L';
% FileQM='test4-2_R';

% FileQM='test4-3_L';
% FileQM='test4-3_R';

%----������ �׽�Ʈ----%
% FileQM='testVT_L';
% FileQM='testVT_R';

%-----------------------------------------------------%
[arrQM, FinalPRNs, FinalTTs] = ReadQM(FileQM);
arrQM = arrQM(find(arrQM(:,3) < 200),:);
% break
yyS = num2str(YY,'%02d');
doyS = num2str(DOY,'%03d');
[gw, gd] = ydoy2gwgd(YY, DOY); %: GPS WEEK ����
[bw, bd] = ydoy2bwbd(YY, DOY); %: BDS WEEK ����

%% ��ǥ ����
% TruePos = [-3041235.578 4053941.677 3859881.013]; % A
% TruePos = [-3041241.741 4053944.143 3859873.640]; % B
% TruePos = [-3041210.419 4053863.515 3859778.262]; % C
% TruePos = [-3041230.128 4053871.023 3859754.445]; % D
% TruePos = [-3012534.7413 ,	4078273.9465 ,	3856561.8977];
% TruePos = [-3012538.7163 4078280.5070 3856568.0046]; %170112 ��õ���� �ʱ���ǥ
TruePos = [-3058799.61420451,4083265.35912516,3814946.87192938]; %170116 ������ �ʱ���ǥ

%% Navigation file load
FileNav = strcat('brdc', num2str(DOY,'%03d'), '0.', num2str(YY,'%02d'), 'n');   %: Navigation RINEX file
eph=ReadEPH(FileNav);

%% ���� �ڵ鸵
LeapSecBDS = 14;

%% brdm���Ͽ��� al/be�� ���� ������ brdc���Ͽ��� ������
al=[0 0 0 0];   be=[0 0 0 0]; 

%% ����� ����ġ ����
g_ObsType = 120; % gps C1
g_ObsType_snr = 141;
c_ObsType = 220; % bds C1
c_ObsType_snr = 241;

%% QM �غ�
QM = SelectQM_gc(arrQM, g_ObsType, c_ObsType);
QM_snr = SelectQM_gc(arrQM, g_ObsType_snr, c_ObsType_snr);
FinalPRNs = unique(QM(:,2));
FinalTTs = unique(QM(:,1));
% break

%% ��ǥ �ʱ�ġ ���� �� ���浵 ��ȯ
AppPos = TruePos;
gd = xyz2gd(AppPos); AppLat = gd(1); AppLon = gd(2);

%% ������ ���� �Ű����� ����
Maxiter = 10;
EpsStop = 1e-4;
ctr = 1; ctr2 = 1;
eleCut = 15;
x = [AppPos ctr]';
%%
NoEpochs = length(FinalTTs);
estm = zeros(NoEpochs, 6);
Scount = zeros(NoEpochs, 1);
nEst = 0;
for j = 1:NoEpochs % 60:65%
    x = [AppPos ctr]';
    for iter = 1:Maxiter
        gs = FinalTTs(j);
        HTH = zeros(4,4);
        HTy = zeros(4,1);
        indexQM = find(QM(:,1) == gs);
        QM1e = QM(indexQM,:);
        QM1e(:,2) = QM1e(:,2) -100;
        QM1e_snr = QM_snr(indexQM,:);
        QM1e_snr(:,2) = QM1e_snr(:,2) -100;
        NoSats = length(QM1e);
        NoSatsUsed = 0;

        if NoSats < 5
            break
        end
        vec_site = x(1:3)';  
        ZHD = TropGPTh(TruePos, gw, gs); %: TROP: GPT
        for i = 1:NoSats
            prn = QM1e(i,2);
            obs = QM1e(i,4);
            snr = QM1e_snr(i,4);
            icol = PickEPH(eph, prn, gs);
            toe = eph(icol, 8); a = eph(icol, 19); b = eph(icol, 20); c = eph(icol, 21); Tgd = eph(icol, 23);
%             STT=0;
            STT = GetSTTbrdc(gs, prn, eph, vec_site); % ��ȣ���޽ð� ���
%             if prn > 100 && prn < 200
                tc = gs - STT ;             % ��ȣ���޽ð� ����
%             elseif prn > 200 && prn < 300
%                 tc = gs - STT- LeapSecBDS;  % ��ȣ���� �ð� ����... bds���� �ݿ�
%             end
            SatPos = GetSatPosNC_GC(eph, icol, tc); % ������ġ ����
            SatPos = RotSatPos(SatPos, STT); % ��������ȿ�� ���
            vec_rho = SatPos - vec_site';
            rho = norm(vec_rho);
            [az,el] = xyz2azel(vec_rho, AppLat, AppLon);
%             if prn > 100 && prn < 200
%             dIono = 0;
%             dTrop = 0;
%             dIono = ionoKlob(al, be, gs, az, el, x(1:3)); % �̿��� ����(Klobuchar ��)
%             dTrop = ZHD2SHD(gw, gs, TruePos, el, ZHD); % ����� ����
%             elseif prn > 200 && prn < 300
%             dIono = 0;
%             dTrop = 0;
            dIono = ionoKlob(al, be, gs, az, el, x(1:3)); % �̿��� ����(Klobuchar ��)
            dTrop = ZHD2SHD(gw, gs, TruePos, el, ZHD); % ����� ����
%             end
%             dRel=0;
            dRel = GetRelBRDC(eph, icol, tc); % ��뼺ȿ��
            dtSat = a + b*(tc - toe) + c*(tc - toe)^2 - Tgd + dRel; % �����ð���� ���... �׷������, ��뼺ȿ�� ����
            if el >=eleCut % �Ӱ����
                W = 1;
%                 W = MakeW_elsnr(el,snr);
%                 if prn > 100 && prn < 200
%                     PRC = PickPRC(RTCM,prn,gs);
                    PRC = 0;
                    com = rho + x(4) - CCC*dtSat + dTrop + dIono - PRC; % gps ��갪
                    H = [-vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 1];
%                 elseif prn > 200
%                     PRC = PickPRC(RTCM,prn,gs);
%                     PRC = 0;
%                     com = rho + x(5) - CCC*dtSat + dTrop + dIono - PRC; % bds ��갪
%                     H = [-vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 0 1];
%                 end
                y = obs - com;
                HTH = HTH + H'*W*H;
                HTy = HTy + H'*W*y;
                NoSatsUsed = NoSatsUsed + 1;
            end
        end
        P = inv(HTH);
        xhat = inv(HTH) * HTy;
        x = x + xhat;
        if norm(xhat) < EpsStop;
            nEst = nEst + 1;
            estm(nEst,1) = gs;
            estm(nEst,2:5) = x(1:4)';
            Scount(nEst,1) = NoSatsUsed;
            fprintf('%8d : %3d : %8.2f : %8.2f\n', j, iter, x(1)' - TruePos(1), x(2)' - TruePos(2));
            break;
        end
    end   
end
estm = estm(1:nEst,:);
Scount=Scount(1:nEst, :);
% figure(101);
[dXYZ, dNEV]=PosErrors(estm(:,1), TruePos,estm(:,2:5),Scount);
break

%% �Է¹��� �ڷ��� ũ�� �� ������ �ش��ϴ� ���浵 ����
NoPos = length(estm(:,1));
gd = xyz2gd(TruePos); TrueLat = gd(1); TrueLon = gd(2); 
%% ������ XYZ�� ���� XYZ�� ���̰� ����
dXYZ = zeros(NoPos,3);
for k = 1:NoPos
    gd(k,:) = xyz2gd(estm(k,2:4));
end
estm_L=estm;
estm_R=estm;
% toc
% break


open ariport_final2.fig
hold on
plot(gd(:,2),gd(:,1),'o:')
axis([126.45219 126.45255 37.44295 37.44325])
xlabel('longitude')
ylabel('latitude')
title('161223 ��õ���� Pole 1')
break

