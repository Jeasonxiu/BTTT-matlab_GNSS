function [estm] = DPP_g(QMfile,eph,TruePos, RTCM, DOY, YY)

%
%   function [estm] = PP_g(QMfile,eph,TruePos, DOY, YY)
%
%   Read the files(QMfile, eph, TruePos, DOY, YY) and estimate Position(GPS/BDS)
%
%   input QMfile
%   input eph matrix : ex> readEPH_all(navfile)
%   input TruePos : 1 X 3 ecef coordinates
%
%   Example : [anything] = PP_gc(QMfile, eph, TruePos, 17, 016)
%
%   coded by Joonseong Gim, Feb 02, 2017
%
%

% clear all;
% close all

% % % %
% QMfile = 'QM170125_A';
% QMfile = 'ublox_joon';
% QMfile = 'ublox_hyunu';
% QMfile = 'QM170322_Bs_1';         % 오송
% QMfile = 'DSDPB_17213_HV';          % 대성 B point
% navfile = 'brdm0250.17p';
% navfile = 'brdm0470.17p';
% navfile = 'brdm0810.17p';         % 오송
% navfile = 'brdm0820.17p';         % 오송
% navfile = 'brdm2130.17p';         % 대성 B point
% % % TruePos = [-3041235.578 4053941.677 3859881.013];
% YY = 17; DOY =025;
% YY = 17; DOY =047;
% YY = 17; DOY =081;        % 오송
% YY = 17; DOY =213;        % 대성 B point
% TruePos = [-3058799.61420451,4083265.35912516,3814946.87192938];
% TruePos = [-3041235.57800000,4053941.67700000,3859881.01300000];    % 대성 A
% TruePos = [-3027386.463213997 4071581.638074351 3852036.292033684]; % 현우집앞
% TruePos = [-3108706.97103816,4078522.84147406,3779757.23816543];    % 오송 1
% TruePos = [-3108697.15866998,4078501.37490046,3779789.12574991];    % 오송 2
% TruePos = [-3041241.741 4053944.143 3859873.640];       % 대성 B point

%% PRC load
% % load('RTCM170802.mat');
load('RTCM170801.mat');
%% 불변 변수 설정: 빛의 속도, 관측치
CCC = 299792458.;   % CCC = Speed of Light [m/s]
%% 임계고도각 설정
eleCut = 15;
%% 윤초 핸들링
LeapSecBDS = 14;

%% 사용할 관측치 선정
g_ObsType = 120; % gps C1
g_ObsType_snr = 141;
c_ObsType = 220; % bds C1
c_ObsType_snr = 241;

%% QM 파일 읽어들여서 행렬로 저장하고, 사용할 관측치 추출
[arrQM, FinalPRNs, FinalTTs] = ReadQM(QMfile);
arrQM = arrQM(find(arrQM(:,3) < 200),:);
arrQM(:,1) = round(arrQM(:,1));

% arrQM_song = arrQM(min(find(arrQM(:,1) == 387531)):max(find(arrQM(:,1) == 387842)),:);
% arrQM_no = arrQM(max(find(arrQM(:,1) == 387842))+1:end,:);
% arrQM_gps = arrQM_song(find(arrQM_song(:,3) < 200),:);
% arrQM_bds = arrQM_song(find(arrQM_song(:,3) > 200),:);
% arrQM_gps = arrQM_gps(find(arrQM_gps(:,2) ~= 1 & arrQM_gps(:,2) ~= 7 & arrQM_gps(:,2) ~= 4 & arrQM_gps(:,2) ~= 30 & arrQM_gps(:,2) ~= 26),:);
% arrQM_bds = arrQM_bds(find(arrQM_bds(:,2) ~= 4),:);
% arrQM_song = [arrQM_gps; arrQM_bds];
% arrQM_song = arrQM_song(find(arrQM_song(:,2) ~= 7 & arrQM_song(:,3) > 200),:);
% arrQM_song = arrQM_song(find(arrQM_song(:,2) ~= 1 & arrQM_song(:,2) ~= 23 &...
%     arrQM_song(:,2) ~= 30 & arrQM_song(:,2) ~= 9),:);

% arrQM_song = arrQM_song(find(arrQM_song(:,2) ~= 30 & arrQM_song(:,2) ~= 3 &...
%     arrQM_song(:,2) ~= 26 & arrQM_song(:,2) ~= 7 & arrQM_song(:,2) ~= 9 & arrQM_song(:,2) ~= 23 & arrQM_song(:,2) ~= 4),:);
% arrQM = [arrQM_song;arrQM_no];
QM = SelectQM_gc(arrQM, g_ObsType, c_ObsType);
QM_snr = SelectQM_gc(arrQM, g_ObsType_snr, c_ObsType_snr);
FinalTTs = unique(QM(:,1));
% FinalTTs = FinalTTs(find(FinalTTs(:,1) > 1000),:);

%% Doppler smoothing을 위한 QM 핸들링
% QMd =qmHandle(arrQM);
% load('eph170216.mat')
% load('eph170322.mat')
% load('eph170801.mat')

% %% eph 생성
% [eph, trashPRN, trashT]=ReadEPH_all(navfile);
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

%% 추정과정 시작
NoEpochs = length(FinalTTs);
estm = zeros(NoEpochs, 6);
Scount = zeros(NoEpochs, 1);
nEst = 0;

tic;
for j = 1:NoEpochs
    % for j = 3000:4000
    gs = FinalTTs(j);
    %     QMe_ = qmHandle(QMd.pickQM(gs,':',':'));
    %     QMe = dopplerSM(QMe_.getQM(),5);
    indexQM = find(QM(:,1) == gs);
    QM1e = QM(indexQM,:);
    QM1e_snr = QM_snr(indexQM,:);
    PRC1e = RTCM(find(RTCM(:,1) == gs),:);
    existprn = intersect(QM1e(:,2), unique(eph(find(eph(:,22) == 0),18)));
    existprn = intersect(existprn,PRC1e(:,2));
    NoSats = length(existprn);
    
    for iter = 1:Maxiter
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
        for i = 1:NoSats
            
            prn = QM1e(find(QM1e(:,2) == existprn(i)),2);
            obs = QM1e(find(QM1e(:,2) == existprn(i)),4);
            snr = QM1e_snr(find(QM1e_snr(:,2) == existprn(i)),4);
            icol = PickEPH(eph, prn, gs);
            toe = eph(icol, 8); a = eph(icol, 19); b = eph(icol, 20); c = eph(icol, 21); Tgd = eph(icol, 23);
            %             STT=0;
            STT = GetSTTbrdm(gs, eph, icol, x(1:3)); % 신호전달시간 계산
            if prn < 300
                if prn > 100 && prn < 200
                    tc = gs - STT ;             % 신호전달시간 보정
                elseif prn > 200 && prn < 300
                    tc = gs - STT- LeapSecBDS;  % 신호전달 시간 보정... bds윤초 반영
                end
                SatPos = GetSatPosNC_GC(eph, icol, tc); % 위성위치 산출
                SatPos = RotSatPos(SatPos, STT); % 지구자전효과 고려
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
            dRel = GetRelBRDC(eph, icol, tc); % 상대성효과
            dtSat = a + b*(tc - toe) + c*(tc - toe)^2 - Tgd + dRel; % 위성시계오차 계산... 그룹딜레이, 상대성효과 보정
            
            if el >=eleCut % 임계고도각
                W = 1;
                %                 W = MakeW_elsnr(el,snr);
                if prn > 100 && prn < 200
                    PRC = PickPRC(PRC1e,prn,gs);
%                     PRC = 0;
%                     com = rho + x(4) - CCC*dtSat + dTrop + dIono - PRC; % gps 계산값
                    com = rho + x(4) - CCC*dtSat - PRC; % gps 계산값
                    H = [-vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 1];
                    NoGPSsUsed = NoGPSsUsed + 1;
                elseif prn > 200 && prn < 300
                    PRC = 0;
                    com = rho + x(5) - CCC*dtSat + dTrop + dIono - PRC; % bds 계산값
                    H = [-vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 1];
                    NoBDSsUsed = NoBDSsUsed + 1;
                end
                y = obs - com;
                HTH = HTH + H'*W*H;
                HTy = HTy + H'*W*y;
                
            end
        end
        if NoGPSsUsed >=4
            P = inv(HTH);
            xhat = inv(HTH) * HTy;
            x = x + xhat;
            if norm(xhat) < EpsStop;
                nEst = nEst + 1;
                estm(nEst,1) = gs;
                estm(nEst,2:5) = x(1:4);
                estm(nEst,7) = NoGPSsUsed;
                estm(nEst,8) = NoBDSsUsed;
                estm(nEst,9) = NoGLOsUsed;
                %             Scount(nEst,1) = NoUsed_gps;
                fprintf('%8d : %3d : %8.2f : %8.2f\n', j, iter, x(1)' - TruePos(1), x(2)' - TruePos(2));
                break;
            end
        else
            break;
        end
    end
end
estm = estm(1:nEst,:);
Scount=Scount(1:nEst, :);
toc;
% [dXYZ, dNEV] = PosTErrors2(estm(:,1), TruePos, estm(:,2:5));
