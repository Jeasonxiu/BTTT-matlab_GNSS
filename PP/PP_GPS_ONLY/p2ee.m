clear all;
close all
%
% % % % %
% % QMfile = 'QM170125_A';
% % QMfile = 'ublox_joon';
% % QMfile = 'ublox_hyunu';
% % QMfile = 'QM170322_Bs_1';         % 오송
% % QMfile = 'DSDPB_17213_HV';          % 대성 B point
% % QMfile = 'ANSG_17214_ud1-1';
QMfile = 'QTEGN_18001';
% % navfile = 'brdm0250.17p';
% % navfile = 'brdm0470.17p';
% % navfile = 'brdm0810.17p';         % 오송
% % navfile = 'brdm0820.17p';         % 오송
% % navfile = 'brdm2130.17p';         % 대성 B point
navfile = 'brdc0010.18n';
% % % % TruePos = [-3041235.578 4053941.677 3859881.013];
% % YY = 17; DOY =025;
% % YY = 17; DOY =047;
% % YY = 17; DOY =081;        % 오송
% % YY = 17; DOY =213;        % 대성 B point
% DOY = 214; YY  = 17;        % 안성W
DOY = 001; YY  = 18;        % 안성W
% % TruePos = [-3058799.61420451,4083265.35912516,3814946.87192938];
% % TruePos = [-3041235.57800000,4053941.67700000,3859881.01300000];    % 대성 A
% % TruePos = [-3027386.463213997 4071581.638074351 3852036.292033684]; % 현우집앞
% % TruePos = [-3108706.97103816,4078522.84147406,3779757.23816543];    % 오송 1
% % TruePos = [-3108697.15866998,4078501.37490046,3779789.12574991];    % 오송 2
% % TruePos = [-3041241.741 4053944.143 3859873.640];       % 대성 B point
% % TruePos = [-3080121.50057624,4058693.37773498,3824124.32021532];    % 안성W 1-1
TruePos=[-3241051.567 4030771.731 3719838.489];
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
QM = SelectQM_gc(arrQM, g_ObsType, c_ObsType);
QM = QM(find(QM(:,2) ~= 102 & QM(:,2) ~=103),:);
% QM_snr = SelectQM_gc(arrQM, g_ObsType_snr, c_ObsType_snr);
FinalTTs = unique(QM(:,1));



%% eph 생성
[eph]=ReadEPH(navfile);
eph(:,18) = eph(:,18) + 100;
% load('eph170322.mat');
% load('eph170801.mat');      % 대성 B point
% load('eph170802.mat');
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
    %     QM1e_snr = QM_snr(indexQM,:);
    existprn = intersect(QM1e(:,2), unique(eph(find(eph(:,22) == 0),18)));
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
        cnt2=1;
        for i = 1:NoSats
            
            prn = QM1e(find(QM1e(:,2) == existprn(i)),2);
            obs = QM1e(find(QM1e(:,2) == existprn(i)),4);
            %             snr = QM1e_snr(find(QM1e_snr(:,2) == existprn(i)),4);
            icol = PickEPH(eph, prn, gs);
            toe = eph(icol, 8); a = eph(icol, 19); b = eph(icol, 20); c = eph(icol, 21); Tgd = eph(icol, 23);
            STT = GetSTTbrdc(eph(icol,:), gs, vec_site);
            
            tc = gs - STT ;             % 신호전달시간 보정
            dRel = GetRelBRDC(eph, tc); % 상대성효과
            dtSat = a + b*(tc - toe) + c*(tc - toe)^2 - Tgd + dRel;
            SatPos = GetSatPosNC_GC(eph, icol, tc); % 위성위치 산출
            SatPos = RotSatPos(SatPos, STT); % 지구자전효과 고려
            vec_rho = SatPos - vec_site;
            rho = norm(vec_rho);
            [az,el] = xyz2azel(vec_rho, AppLat, AppLon);
            
            dIono = klobuchar(al, be, gs, az, el, x(1:3)); % 이온층 보정(Klobuchar 모델)
            dTrop = ZHD2SHD(gw, gs, TruePos, el, ZHD); % 대류권 보정
            %             dIono = 0;
            %             dTrop = 0;
            if el >=eleCut % 임계고도각
                W(cnt2,cnt2) = 1;
                %                 W = MakeW_elsnr(el,snr);
                PRC = 0;
                com = rho + x(4) - CCC*dtSat + dTrop + dIono - PRC; % gps 계산값
                H(cnt2,:) = [-vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 1];
                NoGPSsUsed = NoGPSsUsed + 1;
                y(cnt2,1) = obs - com;
                cnt2=cnt2+1;
                
            end
        end
        if NoGPSsUsed >=4
            HTH = H(1:cnt2-1,:)'*inv(W(1:cnt2-1,1:cnt2-1))*H(1:cnt2-1,:);
            HTy = H(1:cnt2-1,:)'*inv(W(1:cnt2-1,1:cnt2-1))*y(1:cnt2-1,:);
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
                fprintf('%8d : %3d : %8.2f : %8.2f : %8.2f\n', j, iter, x(1)' - TruePos(1), x(2)' - TruePos(2),x(3)' - TruePos(3));
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