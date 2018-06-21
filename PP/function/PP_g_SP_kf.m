% function [estm] = PP_g_SP_kf(QMfile,eph,TruePos, DOY, YY)

%
%   [estm] = PP_g_SP_kf(QMfile,eph,TruePos, DOY, YY)
%
%   Read the files(QMfile, eph, TruePos, DOY, YY) and estimate Position(GPS)
%
%   input QMfile
%   input eph matrix : ex> readEPH_all(navfile)
%   input TruePos : 1 X 3 ecef coordinates
%
%   Example : [anything] = PP_g_SP_kf(QMfile, eph, TruePos, 17, 016)
%
%   coded by Joonseong Gim, Sep 02, 2017
%
%

clear all;
close all

% % % %
% QMfile = 'QM170125_A';
% QMfile = 'ublox_joon';
% QMfile = 'ublox_hyunu';
% QMfile = 'QM170322_Bs_1';         % 오송
% QMfile = 'DSDPB_17213_FH';          % 대성 B point
% QMfile = 'DSDPB_17213_RH';          % 대성 B point
% QMfile = 'DSDPB_17213_HV';          % 대성 B point
% QMfile = 'DSDPB_17213_LV';          % 대성 B point
% QMfile = 'DSDPB_170914_10';          % 대성 B point
% QMfile = 'DSDPB_170919_2';          % 대성 Round
% QMfile = 'DSDPB_171013_1';          % 대성 B point
% navfile = 'brdm0250.17p';
% navfile = 'brdm0470.17p';
% navfile = 'brdm0810.17p';         % 오송
% navfile = 'brdm0820.17p';         % 오송
% navfile = 'brdm2130.17p';   YY = 17; DOY =213;  load('eph170801.mat');      % 대성 B point
% navfile = 'brdm2570.17p';   YY = 17; DOY =257;  load('eph170914.mat');      % 대성 B point Nexus9
navfile = 'brdm2860.17p';   YY = 17; DOY =286;  load('eph171013.mat');      % 대성 B point Nexus9
% navfile = 'brdc2620.17n';   YY = 17; DOY =262;  load('eph170919.mat');      % 대성 B point Nexus9

% TruePos = [-3058799.61420451,4083265.35912516,3814946.87192938];
% TruePos = [-3041235.57800000,4053941.67700000,3859881.01300000];    % 대성 A
% TruePos = [-3027386.463213997 4071581.638074351 3852036.292033684]; % 현우집앞
% TruePos = [-3108706.97103816,4078522.84147406,3779757.23816543];    % 오송 1
% TruePos = [-3108697.15866998,4078501.37490046,3779789.12574991];    % 오송 2
TruePos = [-3041241.741 4053944.143 3859873.640];       % 대성 B point
%% 불변 변수 설정: 빛의 속도, 관측치
CCC = 299792458.;   % CCC = Speed of Light [m/s]
%% 임계고도각 설정
eleCut = 15;
%% 윤초 핸들링
LeapSecBDS = 14;

%% 사용할 관측치 선정
g_ObsType = 120; % gps C1
g_ObsType_snr = 141;
g_ObsType_PrSigmaM = 121;
g_ObsType_prrSigmaMps = 122;
c_ObsType = 220; % bds C1
c_ObsType_snr = 241;

%% QM 파일 읽어들여서 행렬로 저장하고, 사용할 관측치 추출
[arrQM, FinalPRNs, FinalTTs] = ReadQM(QMfile);
arrQM = arrQM(find(arrQM(:,3) < 200),:);
% arrQM(:,1) = round(arrQM(:,1));

%% QM 핸들링
QM = SelectQM_gc(arrQM, g_ObsType, c_ObsType);
QM_snr = SelectQM_gc(arrQM, g_ObsType_snr, c_ObsType_snr);
QM_prr = load(QMfile);
QM_PrSi gmaM = SelectQM_gc(QM_prr, g_ObsType_PrSigmaM, c_ObsType_snr);
FinalTTs = unique(QM(:,1));

% FinalTTs = FinalTTs(find(FinalTTs(:,1) > 1000),:);

%% Doppler smoothing을 위한 QM 핸들링
% QMd =qmHandle(arrQM);
% load('eph170216.mat')
% load('eph170322.mat')

% %% eph 생성
% [eph, trashPRN, trashT]=ReadEPH_all(navfile);
% [eph, trashPRN, trashT]=ReadEPH(navfile);
% load('eph170322.mat');

%% 항법메시지를 읽어들여서 행렬로 저장하고, Klobuchar 모델 추출
gps_nav = strcat('brdc', num2str(DOY,'%03d'), '0.', num2str(YY,'%02d'), 'n');   %: Navigation RINEX file
fid = fopen(gps_nav,'r');
if fid == -1
    al = zeros(4,1); be = zeros(4,1);
else
    [al, be] = GetALBE(gps_nav);
end

%% 항법메시지 파일 이름을 이용해 YY, DOY 생성
[gw, gd] = ydoy2gwgd(YY, DOY); %: GPS WEEK 결정

%% 초기좌표 획득
AppPos = TruePos;
gd = xyz2gd(AppPos); AppLat = gd(1); AppLon = gd(2);

%% 추정에 필요한 초기치 설정
Maxiter = 10;
EpsStop = 1e-5;
ctr = 1; ctr2 = 1;
x = [AppPos ctr]; x = x';
% x = [0, 0, 0, 0]';

%% 추정과정 시작
NoEpochs = length(FinalTTs);
estm = zeros(NoEpochs, 6);
Scount = zeros(NoEpochs, 1);
nEst = 0;

%% Kalman filter
NoPara = 4;         % Kalman filter

% P 초기화 
P = eye(NoPara);            % Kalman filter
P(1:3,1:3) = P(1:3,1:3)*1;  % Kalman filter
% P(4,4) = P(4,4)*0.5;        % Kalman filter
P(4,4) = P(4,4)*0.5;        % Kalman filter

% Q 초기화
Q = eye(NoPara)*0.0003;            % Kalman filter
% Q = eye(NoPara)*0.01;            % Kalman filter
Q(4,4) = 10000;           % Kalman filter

tic;
for j = 1:NoEpochs
    gs = FinalTTs(j);
    if j>2
    gs2 = FinalTTs(j-1);
    indexQM2 =find(QM(:,1) == gs2); 
    QM1et_1 =QM(indexQM2,:);
    gs3 = FinalTTs(j-2);
    indexQM3 =find(QM(:,1) == gs3); 
    QM1et_2 =QM(indexQM2,:);
    end
    indexQM = find(QM(:,1) == gs);
    QM1e = QM(indexQM,:);
    QM1e_snr = QM_snr(indexQM,:);
    QM1e_PrSigmaM = QM_PrSigmaM(indexQM,:);
    existprn = intersect(QM1e(:,2), unique(eph(find(eph(:,22) == 0),18)));
    NoSats = length(existprn);
    
%     for iter = 1:Maxiter
        for iter = 1:1
        HTH = zeros(4,4);
        HTy = zeros(4,1);
        
        if NoSats <= 4
            break
        end
        vec_site = x(1:3)';
        ZHD = TropGPTh(TruePos, gw, gs); %: TROP: GPT
        NoGPSsUsed = 0;
        NoBDSsUsed = 0;
        NoGLOsUsed = 0;
        cnt2=1;
        obs2 = []; obs3=[];
        for i = 1:NoSats
            
            prn = QM1e(find(QM1e(:,2) == existprn(i)),2);
            obs = QM1e(find(QM1e(:,2) == existprn(i)),4);
            if j>2
            obs3 = QM1et_2(find(QM1et_2(:,2) == existprn(i)),4);
            obs2 = QM1et_1(find(QM1et_1(:,2) == existprn(i)),4);
            end
            if ~isempty(obs3) & ~isempty(obs2)
                delPr = ((obs-obs2) + (obs2-obs3))/2;
            elseif ~isempty(obs3)
                delPr = (obs-obs3)/2;
            elseif ~isempty(obs2)
                delPr = (obs-obs2);
            else 
                delPr = 0;
            end
                
                
                
                
            snr = QM1e_snr(find(QM1e_snr(:,2) == existprn(i)),4);
            SNR(cnt2,cnt2) = snr;
            icol = PickEPH(eph, prn, gs);
            toe = eph(icol, 8); a = eph(icol, 19); b = eph(icol, 20); c = eph(icol, 21); Tgd = eph(icol, 23); Toc = eph(icol, 26);
            STT = obs/CCC;
            
            if prn < 300
                if prn > 100 && prn < 200
                    tc = gs - STT ;             % 신호전달시간 보정
                elseif prn > 200 && prn < 300
                    tc = gs - STT- LeapSecBDS;  % 신호전달 시간 보정... bds윤초 반영
                end
                dRel = GetRelBRDC(eph, icol, tc); % 상대성효과
                dtSat = a + b*(tc - Toc) + c*(tc - Toc)^2 - Tgd + dRel; % 위성시계오차 계산... 그룹딜레이, 상대성효과 보정
                SatPos = GetSatPosNC_GC(eph, icol, tc-dtSat); % 위성위치 산출
                SatPos = RotSatPos(SatPos, STT-dtSat); % 지구자전효과 고려
            end
            vec_rho = SatPos - vec_site';
            rho = norm(vec_rho);
            [az,el] = xyz2azel(vec_rho, AppLat, AppLon);
            
            if prn > 100 && prn < 200
                dIono = klobuchar(al, be, gs, az, el, x(1:3)); % 이온층 보정(Klobuchar 모델)
                dTrop = ZHD2SHD(gw, gs, TruePos, el, ZHD); % 대류권 보정
            elseif prn > 200 && prn < 300
                dIono = ionoKlob(al, be, gs, az, el, x(1:3)); % 이온층 보정(Klobuchar 모델)
                dTrop = ZHD2SHD(gw, gs, TruePos, el, ZHD); % 대류권 보정
            end
%             dIono = 0;
%             dTrop = 0;
            %             if el >=eleCut % 임계고도각
            if el >=0 % 임계고도각
                W(cnt2,cnt2) = 1;
%                 W(cnt2,cnt2) = 1/QM1e_PrSigmaM(find(QM1e_PrSigmaM(:,2) == prn),4);
% if j >2
%                 W(cnt2,cnt2) = abs(delPr)^2;
% end
                W(cnt2,cnt2) = QM1e_PrSigmaM(find(QM1e_PrSigmaM(:,2) == prn),4);
                
                
                if prn > 100 && prn < 200
                    PRC = 0;
                    com = rho + x(4) - CCC*dtSat + dTrop + dIono - PRC; % gps 계산값
                    H(cnt2,:) = [-vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 1];
                    NoGPSsUsed = NoGPSsUsed + 1;
                elseif prn > 200 && prn < 300
                    PRC = 0;
                    com = rho + x(5) - CCC*dtSat + dTrop + dIono - PRC; % bds 계산값
                    H(cnt2,:) = [-vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 1];
                    NoBDSsUsed = NoBDSsUsed + 1;
                end
                y(cnt2,1) = obs - com;
                cnt2=cnt2+1;
            end
        end
        if NoGPSsUsed >=4
            xp = x;             % Kalman filter
            Po = P;             % Kalman filter
            Pp = P + Q;         % Kalman filter
%             K = Pp*H(1:cnt2-1,:)'/(H(1:cnt2-1,:)*Pp*H(1:cnt2-1,:)'+W(1:cnt2-1,1:cnt2-1)*100);    % Kalman filter
            K = Pp*H(1:cnt2-1,:)'*inv(H(1:cnt2-1,:)*Pp*H(1:cnt2-1,:)'+W(1:cnt2-1,1:cnt2-1));    % Kalman filter
            x = xp + K*y(1:cnt2-1,:);       % Kalman filter
            P = Pp - K*H(1:cnt2-1,:)*Pp;    % Kalman filter
%             HTH = pinv(W(1:cnt2-1,1:cnt2-1)*H(1:cnt2-1,:));
%             xhat = HTH * W(1:cnt2-1,1:cnt2-1)*y(1:cnt2-1,:);
%             x = x + xhat;
%             if norm(xhat) < EpsStop;
                nEst = nEst + 1;
                estm(nEst,1) = round(gs);
                estm(nEst,2:5) = x(1:4);
                estm(nEst,7) = NoGPSsUsed;
                estm(nEst,8) = NoBDSsUsed;
                estm(nEst,9) = NoGLOsUsed;
                %             Scount(nEst,1) = NoUsed_gps;
                estm_gd(nEst,:) = xyz2gd(x(1:3));
                fprintf('%8d : %3d : %8.2f : %8.2f : %8.2f\n', j, iter, x(1)' - TruePos(1), x(2)' - TruePos(2),x(3)' - TruePos(3));
%                 break;
                
%             end
%         else
%             break;
        end
    end
end
estm = estm(1:nEst,:);
Scount=Scount(1:nEst, :);
toc;


load([QMfile,'_gg.mat']);
load([QMfile,'_fix.mat']);
estm_fix(:,1) = round(estm_fix(:,1) -32400);
estm_fix = [estm_fix, zeros(length(estm_fix),5)];
estm_fix(:,5:7) = estm_gg(1:length(estm_fix),5:7);
estm_gg=estm_fix;
% estm_gg(:,1:4) = estm_fix(1:length(estm_gg(:,1)),1:4);
% estm_gg(:,1) = estm(:,1);

[dXYZ, dNEV] = PosTErrorsSP(estm(:,1), TruePos, estm(:,2:5),[estm(:,1),estm(:,7)]);
[dXYZ_gg, dNEV_gg] = PosTErrorsSP(estm_gg(:,1), TruePos, estm_gg(:,2:5),[estm_gg(:,1),estm_gg(:,7)]);
% estm_GG = estm_gg(301:end,:);
% ESTM = estm(find(estm(:,1) == min(estm_GG(:,1))):find(estm(:,1) == max(estm_GG(:,1))),:);
% [dXYZ, dNEV] = PosTErrorsSP(ESTM(:,1), TruePos, ESTM(:,2:5),[ESTM(:,1),ESTM(:,7)]);
% [dXYZ_gg, dNEV_gg] = PosTErrorsSP(estm_GG(:,1), TruePos, estm_GG(:,2:5),[estm_GG(:,1),estm_GG(:,7)]);
% estm_gd = xyz2gd(median(estm(:,2:4)))
% Tru_gd = xyz2gd(TruePos)
% figure(3)
% plot(estm_gd(:,2),estm_gd(:,1),'bo');
% plot_google_map