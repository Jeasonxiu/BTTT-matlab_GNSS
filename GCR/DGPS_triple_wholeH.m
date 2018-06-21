warning off;
clear all; close all;

%% 자료처리 원시데이터
% % QMfile = 'QM170125_A';
% % QMfile = 'QM170407m';
% QMfile = 'QM170616b';
% % navfile = 'brdm0250.17p';
% % navfile = 'brdm0970.17p';
% navfile = 'brdm1670.17p';
% 
% %% PRC
% % PRC_all = load('PPS1_170407.t41');
% PRC_all = load('PPS1_170616.t41');
% %% 기본 파라미터
% % TruePos = [-3041235.578 4053941.677 3859881.013];       % 대성 A point
% % TruePos = [-3003051.386, 4059906.022, 3883100.080];        % 강화 기준국
% TruePos = [-3041241.74100000,4053944.14300000,3859873.64000000];        % 대성 B
% % YY = 17; DOY =025;
% % YY = 17; DOY =097;
% YY = 17; DOY =167;
% 
% %% NMEA
% % NMEAfile = 'PPS1_170407.NMEA';
% NMEAfile = '170616_2.cnb';
% % % [nmea] = writeNMEA(NMEAfile,2017,04,07);
% [nmea] = writeNMEA(NMEAfile,2017,06,16);
% 
% %% 불변 변수 설정: 빛의 속도, 관측치
% CCC = 299792458.;   % CCC = Speed of Light [m/s]
% %% 임계고도각 설정
% eleCut = 15;
% %% 윤초 핸들링
% LeapSecBDS = 14;        % BDS
% LeapSec = 18;           % GLO
% 
% %% 사용할 관측치 선정
% g_ObsType = 120; % gps C1
% c_ObsType = 220; % bds C1
% r_ObsType = 320; % bds C1
% 
% %% QM 파일 읽어들여서 행렬로 저장하고, 사용할 관측치 추출
% [arrQM, FinalPRNs, FinalTTs] = ReadQM(QMfile);
% QM = SelectQM_gcr(arrQM, g_ObsType, c_ObsType, r_ObsType);
% if isempty(find(unique(arrQM(:,3)) == 141))
%     QMsnr = [];
%     nosnr = 1;
% else
%     QMsnr = SelectQM_gcr(arrQM, 141, 241, 341);
%     nosnr = 0;
% end
% FinalPRNs = unique(QM(:,2));
% FinalTTs = unique(QM(:,1));
% 
% %% Navigation file load
% FileNav = strcat('brdm', num2str(DOY,'%03d'), '0.', num2str(YY,'%02d'), 'p');   %: Navigation RINEX file
% %% 항법메시지 파일 이름을 이용해 YY, DOY 생성
% [gw, gd] = ydoy2gwgd(YY, DOY); %: GPS WEEK 결정
% TauC = ReadTauC2(FileNav);
% 
% %% eph 생성
% [eph, trashPRN, trashT]=ReadEPH_all(FileNav);
% ephGLO = eph;
% ephGLO=ephGLO(ephGLO(:,1)<400&ephGLO(:,1)>300,:);ephGLO(:,1)=ephGLO(:,1);
% 
% %% 항법메시지를 읽어들여서 행렬로 저장하고, Klobuchar 모델 추출
% gps_nav = strcat('brdc', num2str(DOY,'%03d'), '0.', num2str(YY,'%02d'), 'n');   %: Navigation RINEX file
% fid = fopen(gps_nav,'r');
% if fid == -1
%     al = zeros(4,1); be = zeros(4,1);
% else
%     [al, be] = GetALBE(gps_nav);
% end
% 
% %% 좌표 초기치 설정 및 위경도 변환
% AppPos = TruePos;
% gd = xyz2gd(AppPos); AppLat = gd(1); AppLon = gd(2);
% %% 추정을 위한 매개변수 설정
% Maxiter = 5;
% EpsStop = 1e-4;
% ctr = 1; ctr2 = 1; ctr3 = 1;
% eleCut = 15; deltat = 60;
% x = [AppPos ctr ctr2 ctr3]';          % All system
% 
% %% 연산 파라미터
% NoEpochs = length(FinalTTs);
% SatPosArr_before=zeros(24,25);icolArr_before=zeros(24,2);
% nEst = 0;
% 
% % load('DGPS_triple_170407m.mat');
load('DGPS_triple_170616b.mat');
AppPos = [-3041229.6208  4053927.3693  3859858.8412];
%% 1=GPS, 2=GLO, 3=BDS, 4= GPS/GLO, 5=GPS/BDS, 6=BDS/GLO, 7= triple
select = 2;
nosnr = 1;
switch select
    case 1
        QM = QM(find(QM(:,3) < 200),:);
        QMsnr = QMsnr(find(QMsnr(:,3) < 200),:);
        x = x(1:4,1);
    case 2
        QM = QM(find(QM(:,3) > 300),:);
        QMsnr = QMsnr(find(QMsnr(:,3) > 300),:);
        deltat = 30;
        x = x(1:4,1);
        eleCut = 5;
    case 3
        QM = QM(find(QM(:,3) > 200 & QM(:,3) < 300),:);
        QMsnr = QMsnr(find(QMsnr(:,3) > 200 & QMsnr(:,3) < 300),:);
        x = x(1:4,1);
    case 4
        QM = QM(find(QM(:,3) < 200 | QM(:,3) > 300),:);
        QMsnr = QMsnr(find(QMsnr(:,3) < 200 | QMsnr(:,3) > 300),:);
        x = x(1:5,1);
    case 5
        QM = QM(find(QM(:,3) < 300),:);
        QMsnr = QMsnr(find(QMsnr(:,3) < 300),:);
        x = x(1:5,1);
    case 6
        QM = QM(find(QM(:,3) > 200),:);
        QMsnr = QMsnr(find(QMsnr(:,3) > 200),:);
        x = x(1:5,1);
end

tic;
for j = 1:NoEpochs
    
    gs = FinalTTs(j);
    indexQM = find(QM(:,1) == gs);
    QM1e = QM(indexQM,:);
    if nosnr == 0
        indexQM = find(QMsnr(:,1) == gs);
        QM1e_snr = QMsnr(indexQM,:);
    end
    PRC_1e = PRC_all(find(PRC_all(:,1) == gs),:);
    existprn = intersect(QM1e(:,2), unique(eph(find(eph(:,22) == 0),18)));
    existprn = intersect(existprn, PRC_1e(:,2));
    NoSats = length(existprn);
    
    for iter = 1:Maxiter
        cnt1=1;
        NoSatsUsed = 0;
        NoGPSsUsed = 0;
        NoBDSsUsed = 0;
        NoGLOsUsed = 0;
        cnt2=1;
        if NoSats <= 5
            break
        end
        vec_site = x(1:3)';
        ZHD = TropGPTh(TruePos, gw, gs); %: TROP: GPT
        for i = 1:NoSats
            
            prn = QM1e(find(QM1e(:,2) == existprn(i)),2);
            obs = QM1e(find(QM1e(:,2) == existprn(i)),4);

            if nosnr == 0
                snr = QM1e_snr(find(QM1e_snr(:,2) == existprn(i)),4); end
            if prn < 300
                icol = PickEPH(eph, prn, gs);
                toe = eph(icol, 8); a = eph(icol, 19); b = eph(icol, 20); c = eph(icol, 21); Tgd = eph(icol, 23);
                STT = GetSTTbrdm(gs, eph, icol, x(1:3)); % 신호전달시간 계산
                
                if prn > 100 && prn < 200
                    tc = gs - STT ;             % 신호전달시간 보정
                elseif prn > 200 && prn < 300
                    tc = gs - STT- LeapSecBDS;  % 신호전달 시간 보정... bds윤초 반영
                end
                
                SatPos = GetSatPosNC_GC(eph, icol, tc); % 위성위치 산출
                SatPos = RotSatPos(SatPos, STT); % 지구자전효과 고려
                
            elseif prn > 300
                icol = PickEPH_GLO2(ephGLO, prn, gs);
                icolArr(i,1)=prn; icolArr(i,2)=icol;
                TauN = ephGLO(icol,12); % : tau & gamma - 시계오차 보정에 필요
                GammaN = ephGLO(icol,13);
                ch_num = ephGLO(icol,16); % : channel number - 전리층 보정에 필요
%                 STT = GetSTTbrdc_GLO_ver3(gs, ephGLO, icol, x(1:3), deltat); % 이전 epoch값으로 신호전달시간 계산
                STT = GetSTTbrdc_GLO(gs, ephGLO, icol, x(1:3), deltat); % 이전 epoch값으로 신호전달시간 계산
                tc = gs - STT;
                tc = tc - LeapSec - TauC;
%                 [SatPos,SatVel,SatLS] = GetSatPosGLO_ver3(ephGLO,icol,tc,deltat); % 이전 epoch값으로 위성위치 계산
                [SatPos,SatVel,SatLS] = GetSatPosGLO(ephGLO,icol,tc,deltat); % 이전 epoch값으로 위성위치 계산
                SatPos = RotSatPos(SatPos, STT); % 지구자전효과 고려
            end
            vec_rho = SatPos - vec_site';
            rho = norm(vec_rho);
            [az,el] = xyz2azel(vec_rho, AppLat, AppLon);
            
            if prn > 100 && prn < 200
                dIono = klobuchar(al, be, gs, az, el, x(1:3)); % 이온층 보정(Klobuchar 모델)
                dTrop = ZHD2SHD(gw, gs, TruePos, el, ZHD); % 대류권 보정
                dRel = GetRelBRDC(eph, icol, tc); % 상대성효과
                dtSat = a + b*(tc - toe) + c*(tc - toe)^2 - Tgd + dRel; % 위성시계오차 계산... 그룹딜레이, 상대성효과 보정
            elseif prn > 200 && prn < 300
                dIono = klobuchar(al, be, gs, az, el, x(1:3)); % 이온층 보정(Klobuchar 모델)
                dTrop = ZHD2SHD(gw, gs, TruePos, el, ZHD); % 대류권 보정
                dRel = GetRelBRDC(eph, icol, tc); % 상대성효과
                dtSat = a + b*(tc - toe) + c*(tc - toe)^2 - Tgd + dRel; % 위성시계오차 계산... 그룹딜레이, 상대성효과 보정
            elseif prn > 300 && prn < 400
                dIono = klobuchar(al, be, gs, az, el, x(1:3)); % 이온층 보정(Klobuchar 모델)
                dTrop = ZHD2SHD(gw, gs, TruePos, el, ZHD); % 대류권 보정
                dRel = (-2/CCC^2) * dot(SatPos, SatVel); % 상대론적 효과
                tb = ephGLO(icol,2) + LeapSec;
                dtSat = TauN - GammaN*(tc-tb) + TauC + dRel;
            end
                        
            if el >=eleCut % 임계고도각
                if iter == 1
                    PRNS(j,1) = gs;
                    PRNS(j,cnt1+1) = prn;
                    cnt1=cnt1+1;
                end
                if nosnr == 0;
                    elsnr(cnt2,1:3) = [prn, el, snr]; 
                else
                    el(cnt2,1:2) = [prn, el];
                    W(cnt2,cnt2) = MakeW_elpr(el(cnt2,2));
                end
                if prn > 100 && prn < 200
                    PRC = 0;
                    PRC = PRC_1e(find(PRC_1e(:,2) == prn),3);
%                     com = rho + x(4) - CCC*dtSat + dTrop + dIono - PRC; % gps 계산값
                    com = rho + x(4) - CCC*dtSat - PRC; % DGPS 계산값
                    if select == 1
                        H(cnt2,:) = [-vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 1];
                    elseif select == 3 | select == 4 | select == 5
                        H(cnt2,:) = [-vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 1 0];
                    elseif select == 7
                        H(cnt2,:) = [-vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 1 0 0];    % triple
                    end
                    NoGPSsUsed = NoGPSsUsed + 1;
                elseif prn > 200 && prn < 300
                    PRC = 0;
                    PRC = PRC_1e(find(PRC_1e(:,2) == prn),3);
                    if select == 3
%                         com = rho + x(4) - CCC*dtSat + dTrop + dIono - PRC; % bds 계산값
                        com = rho + x(4) - CCC*dtSat - PRC; % DBDS 계산값
                        H(cnt2,:) = [-vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 1];
                    elseif select == 5
%                         com = rho + x(5) - CCC*dtSat + dTrop + dIono - PRC; % bds 계산값
                        com = rho + x(5) - CCC*dtSat - PRC; % DBDS 계산값
                        H(cnt2,:) = [-vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 0 1];
                    elseif select == 6
%                         com = rho + x(4) - CCC*dtSat + dTrop + dIono - PRC; % bds 계산값
                        com = rho + x(4) - CCC*dtSat - PRC; % DBDS 계산값
                        H(cnt2,:) = [-vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 1 0];
                    elseif select == 7
%                         com = rho + x(5) - CCC*dtSat + dTrop + dIono - PRC; % bds 계산값
                        com = rho + x(5) - CCC*dtSat - PRC; % DBDS 계산값
                        H(cnt2,:) = [-vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 0 1 0];    % triple
                    end
                    NoBDSsUsed = NoBDSsUsed + 1;
                elseif prn > 300 && prn < 400
                    PRC = 0;
                    PRC = PRC_1e(find(PRC_1e(:,2) == prn),3);
                    if select == 2
%                         com = rho + x(4) - CCC*dtSat + dTrop + dIono - PRC;
                        com = rho + x(4) - CCC*dtSat - PRC;     % DGLO
                        H(cnt2,:) = [-vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 1];
                    elseif select == 4
%                         com = rho + x(5) - CCC*dtSat + dTrop + dIono - PRC;
                        com = rho + x(5) - CCC*dtSat - PRC;     %DGLO
                        H(cnt2,:) = [-vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 0 1];
                    elseif select == 6
%                         com = rho + x(5) - CCC*dtSat + dTrop + dIono -PRC;
                        com = rho + x(5) - CCC*dtSat - PRC;      % DGLO
                        H(cnt2,:) = [-vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 0 1];
                    elseif select == 7
%                         com = rho + x(6) - CCC*dtSat + dTrop + dIono -PRC;
                        com = rho + x(6) - CCC*dtSat - PRC;     % DGLO
                        H(cnt2,:) = [-vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 0 0 1];    % triple
                    end
                    NoGLOsUsed = NoGLOsUsed + 1;
                end
                y = obs - com;
                Y(cnt2,1) = y;
                NoSatsUsed = NoSatsUsed + 1;
                cnt2= cnt2 + 1;
                if nosnr == 2
                    W(cnt2,cnt2) = MakeW_elpr(el(2));end
            end
        end
        if NoSatsUsed >= 5
            
            %% weighting
            if nosnr == 0
                W = MakeW_elsnr_v1(elsnr(:,2),elsnr(:,3),4); 
            elseif nosnr == 1
                W = eye(cnt2-1);
            end
            
            HTH = H(1:cnt2-1,:)'*inv(W(1:cnt2-1,1:cnt2-1))*H(1:cnt2-1,:);
            HTy = H(1:cnt2-1,:)'*inv(W(1:cnt2-1,1:cnt2-1))*Y(1:cnt2-1,:);
            xhat = inv(HTH) * HTy;
            x = x + xhat;
            if norm(xhat) < EpsStop;
                nEst = nEst + 1;
                estm(nEst,1) = gs;
                if select == 1 | select == 2 | select == 3
                    estm(nEst,2:5) = x(1:4);
                elseif select == 4 | select == 5 | select == 6
                    estm(nEst,2:6) = x(1:5);
                elseif select == 7
                    estm(nEst,2:7) = x(1:6);
                end
                
                estm(nEst,8) = NoGPSsUsed;
                estm(nEst,9) = NoBDSsUsed;
                estm(nEst,10) = NoGLOsUsed;
                estm(nEst,11:12) = [NoSats, NoSatsUsed];      % 관측 위성 수와 사용 위성
                fprintf('%8d : %3d : %8.2f : %8.2f : %8.2f\n', j, iter, x(1)' - TruePos(1), x(2)' - TruePos(2), x(3)' - TruePos(3));
                break;
            end
        else
            break;
        end
    end
end

%% 측위오차 분석 & 그래프 작성
estm = estm(1:nEst, :);
toc;
[dXYZ, dNEV] = PosTErrors2(estm(:,1), TruePos, estm(:,2:5),estm(:,11:12));

subplot(3,4,[1,2,5,6])
hold on; grid on;
plot([-20, 20],[0,0],'-r');
plot([0, 0],[-20,20],'-r');

for i=1:length(estm(:,1))
    estmgd(i,:) = xyz2gd(estm(i,2:4));
end
figure(99)
hold on; grid on;
plot(nmea(:,3), nmea(:,2),'bo')
plot(estmgd(:,2), estmgd(:,1),'ro')
if DOY == 167
    plot(126.876985221076, 37.4797132062492, 'g*','Markersize', 10); end;
if DOY == 167
legend('NMEA','estm','B point');
else legend('NMEA','estm'); end
plot_google_map