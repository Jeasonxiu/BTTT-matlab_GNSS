%% 코드의사거리 GPS/GLO 통합 이중차분 알고리즘
% 31/08/2016 : Joonseong
close all; clear all;
tic
%
%% 불변 변수 설정 : 빛의 속도, 관측치
CCC = 299792458.;   % CCC = Speed of Light [m/s]
GPSObsTypeC1 = 120;      % 사용할 관측지 설정 - 120 : C/A = C1
GPSObsTypeSNR = 141;      % 사용할 관측치 설정 - 141: snr = S1
GLOObsTypeC1 = 320;      % 사용할 관측지 설정 - 120 : C/A = C1
GLOObsTypeSNR = 341;      % 사용할 관측치 설정 - 141: snr = S1

%% System 선택 : 1 = GPS, 2 = GLO, 3 = GPS/GLO
System = 2;

%% 임계고도각 설정
eleCut = 15;

%% True Distance
Truedis = 1.41;         % 신항대로 앞뒤
% Truedis = 0.8;         % 신항대로 좌우

%% QMfile load
% bsQMfile = 'QMjfr5';        % 신항대로 앞뒤
% rvQMfile = 'QMjrr5';        % 신항대로 앞뒤
bsQMfile = 'QJF04_ob219';        % 자바드 신항대로 앞뒤
rvQMfile = 'QJR04_ob219';        % 자바드 신항대로 앞뒤
% bsQMfile = 'QUF04_ob219';        % 유블록 신항대로 좌우
% rvQMfile = 'QUR04_ob219';        % 유블록 신항대로 좌우

%% Base 좌표 계산 및 초기좌표 추출 후 위경도 변환
obsfileBs = 'jf042190.obs';                                 % 자바드 Reference Obsevation
navfileBs = 'jf042190.nav';                                 % 자바드 Reference Navigation
% obsfileBs = 'uf042190g.obs';                                 % 유블록 Reference Obsevation
% navfileBs = 'uf042190g.nav';                                 % 유블록 Reference Navigation
% Bs = PPwC(obsfileBs,navfileBs);                                % without Correction
Bs = PP(obsfileBs,navfileBs);                                % without Correction
AppPos = GetAppPos(obsfileBs);
gd = xyz2gd(AppPos); AppLat = gd(1); AppLon = gd(2);

%% 기준위성 및 제외 PRN 설정
RS = 23;
OPRN = 0;      % 제외 PRN 없는 경우 0

%% QM 파일 읽어들여서 행렬로 저장하고, GPS/GLO PRN 중복 방지를 위해 GLO PRN에 100을 더함
[arrQM1, FinalPRNs1, FinalTTs1] = ReadQM(bsQMfile);         % Base
[arrQM2, FinalPRNs2, FinalTTs2] = ReadQM(rvQMfile);         % Rover
% idxbsglo = find(arrQM1(:,3) > 300); arrQM1(idxbsglo, 2) = arrQM1(idxbsglo, 2) + 100;
% idxrvglo = find(arrQM2(:,3) > 300); arrQM2(idxrvglo, 2) = arrQM2(idxrvglo, 2) + 100;

%% arrQM에서 type 별 QM matrix 생성
QMGPSC1bs = SelectQM(arrQM1, GPSObsTypeC1); QMGPSC1bs = QMGPSC1bs(find(QMGPSC1bs(:,2) ~= OPRN),:);
QMGPSSNRbs = SelectQM(arrQM1, GPSObsTypeSNR); QMGPSSNRbs = QMGPSSNRbs(find(QMGPSSNRbs(:,2) ~= OPRN),:);
QMGPSDOPbs = SelectQM(arrQM1, 131); QMGPSDOPbs = QMGPSDOPbs(find(QMGPSDOPbs(:,2) ~= OPRN),:);
QMGLOC1bs = SelectQM(arrQM1, GLOObsTypeC1);
QMGLOSNRbs = SelectQM(arrQM1, GLOObsTypeSNR); QMGLOSNRbs = QMGLOSNRbs(find(QMGLOSNRbs(:,2) ~= OPRN),:);
QMGLODOPbs = SelectQM(arrQM1, 331); QMGLODOPbs = QMGLODOPbs(find(QMGLODOPbs(:,2) ~= OPRN),:);
QMGPSC1rv = SelectQM(arrQM2, GPSObsTypeC1); QMGPSC1rv = QMGPSC1rv(find(QMGPSC1rv(:,2) ~= OPRN),:);
QMGPSSNRrv = SelectQM(arrQM2, GPSObsTypeSNR); QMGPSSNRrv = QMGPSSNRrv(find(QMGPSSNRrv(:,2) ~= OPRN),:);
QMGPSDOPrv = SelectQM(arrQM2, 131); QMGPSDOPrv = QMGPSDOPrv(find(QMGPSDOPrv(:,2) ~= OPRN),:);
QMGLOC1rv = SelectQM(arrQM2, GLOObsTypeC1); QMGLOC1rv = QMGLOC1rv(find(QMGLOC1rv(:,2) ~= OPRN),:);
QMGLOSNRrv = SelectQM(arrQM2, GLOObsTypeSNR); QMGLOSNRrv = QMGLOSNRrv(find(QMGLOSNRrv(:,2) ~= OPRN),:);
QMGLODOPrv = SelectQM(arrQM2, 331); QMGLODOPrv = QMGLODOPrv(find(QMGLODOPrv(:,2) ~= OPRN),:);

%% GPS/GLO 통합 QM matrix 생성
FinalQMbs = [QMGPSC1bs; QMGLOC1bs];                     % Base GPS/GLO C1
FinalQM_snrbs = [QMGPSSNRbs; QMGLOSNRbs];               % Base GPS/GLO snr
FinalQM_dopbs = [QMGPSDOPbs; QMGLODOPbs];               % Base GPS/GLO snr

FinalQMrv = [QMGPSC1rv; QMGLOC1rv];                     % Rover GPS/GLO C1
FinalQM_snrrv = [QMGPSSNRrv; QMGLOSNRrv];               % Rover GPS/GLO snr
FinalQM_doprv = [QMGPSDOPrv; QMGLODOPrv];               % Rover GPS/GLO snr

%% PRN 제외할 구역 설정
OSTART = 0;
OSTOP = 560445;

%% YY, DOY, GPS week, GPS week day 계산
[YY, DOY] = obs2YYDOY(obsfileBs);
[gw, gwd] = ydoy2gwgd(YY, DOY);

%% YY, DOY 에 맞는 BRDC 파일 로드
FileNavGPS = strcat('brdc',num2str(DOY),'0.',num2str(YY),'n');
FileNavGLO = strcat('brdc',num2str(DOY),'0.',num2str(YY),'g');

%% 항법메시지를 읽어들여서 matrix 생성, Klobuchar 모델, TauC, LeapSec 추출
EphGPS = ReadEPH(FileNavGPS);
[al, be] = GetALBE(FileNavGPS);
EphGLO = ReadEPH_GLO(FileNavGLO);
TauC = ReadTauC(FileNavGLO); %%
LeapSec = GetLeapSec(FileNavGLO); %%

%% 두 QM 파일에서 공통시각(epoch) 추출
FinalTTs = intersect(FinalQMbs(:, 1), FinalQMrv(:, 1));
FinalTTs = intersect(Bs(:, 1), FinalTTs(:, 1));

%% 추정에 필요한 초기치 설정
MaxIter = 10;           % Iterarion 횟수
EpsStop = 1e-5;         % 계산 임계값
ctr = 1; deltat = 1;
x = [AppPos]; x = x';       % 계산 좌표 초기값 Matrix

%% 추정과정 시작
NoEpochs = length(FinalTTs);
EstPos = zeros(NoEpochs,5);
nEst = 0;
cntgps = 0; cntglo = 0;
j=1;
estm = zeros(NoEpochs,6);           % Estimation Matrix
visiSat = zeros(NoEpochs,2);        % GPS/GLO visible Satellites
GPS0_c = zeros(NoEpochs,32);        % GPS O-C matrix
GLOo_c = zeros(NoEpochs,24);        % GLO O-C matrix
% load('DDgpsglotest.mat');

switch System
    case 1
                for j = 1:NoEpochs
            % for j = 1:1
            gs = FinalTTs(j);
            idxQMbs = find(FinalQMbs(:,1) == gs);
            idxQMrv = find(FinalQMrv(:,1) == gs);
            QM_c1bs = FinalQMbs(idxQMbs,:);             % Base GPS/GLO C1, 1 Epoch
            QM_snrbs = FinalQM_snrbs(idxQMbs,:);        % Base GPS/GLO SNR, 1 Epoch
            QM_dopbs = FinalQM_dopbs(idxQMbs,:);        % Base GPS/GLO DOP, 1 Epoch
            QM_c1rv = FinalQMrv(idxQMrv,:);             % Rover GPS/GLO C1, 1 Epoch
            QM_snrrv = FinalQM_snrrv(idxQMrv,:);        % Rover GPS/GLO SNR, 1 Epoch
            QM_doprv = FinalQM_doprv(idxQMrv,:);        % Rover GPS/GLO DOP, 1 Epoch
            
            GpsQMbs = QM_c1bs(find(QM_c1bs(:,3) == 120),:); GpsSatsbs = length(GpsQMbs(:,2));               % Base GPS C1
            GloQMbs = QM_c1bs(find(QM_c1bs(:,3) == 320),:); GloSatsbs = length(GloQMbs(:,2));               % Base GLO C1
            GpsQM_snrbs = QM_snrbs(find(QM_snrbs(:,3) == 141),:); GpsSatsbs = length(GpsQM_snrbs(:,4));     % Base GPS SNR
            GloQM_snrbs = QM_snrbs(find(QM_snrbs(:,3) == 341),:); GloSatsbs = length(GloQM_snrbs(:,4));     % Base GLO SNR
            
            GpsQMrv = QM_c1rv(find(QM_c1rv(:,3) == 120),:); GpsSatsrv = length(GpsQMrv(:,2));               % Rover GPS C1
            GloQMrv = QM_c1rv(find(QM_c1rv(:,3) == 320),:); GloSatsrv = length(GloQMrv(:,2));               % Rover GLO C1
            GpsQM_snrrv = QM_snrrv(find(QM_snrrv(:,3) == 141),:); GpsSatsrv = length(GpsQM_snrrv(:,4));     % Rover GPS SNR
            GloQM_snrrv = QM_snrrv(find(QM_snrrv(:,3) == 341),:); GloSatsrv = length(GloQM_snrrv(:,4));     % Rover GLO SNR
            
            %% 해당 시각 Bs의 위치 값 찾기
            Base(j,:) = Bs(find(Bs(:,1) == gs),:);
            TruePosBs_PP = Base(j,2:4);
%             TruePosBs_PP = [-3041241.741 4053944.143 3859873.640];
            base_gd(j,1) = gs;
            base_gd(j,2:4) = xyz2gd(TruePosBs_PP(1:3)); % Base의 xyz값 gd 로 변환
            AppLatBs = base_gd(j,2); AppLonBs = base_gd(j,3);
            vec_site = x(1:3)';
            
            %% GLO 위성 좌표 계산을 위한 적분 초기값 생성
            sT= mod(gs,86400)/3600;
            sTh = floor(sT); sTm = sT - sTh;
            if j == 1 % 시작시간 위성 위치 array
                %         [Sat_ar] = GetSatPosGLO_new(EphGlo,gs,deltat);
                [Sat_ar] = GetSatPosGLO_my(EphGLO,gs,deltat);
                %         fprintf('시작 시간 위치 계산\n');
            elseif (j ~= 1 && sTm == 0) % 정각일 때
                clear Sat_ar
                %         [Sat_ar] = GetSatPosGLO_new(EphGlo,gs,deltat);
                [Sat_ar] = GetSatPosGLO_my(EphGLO,gs,deltat);
                %         fprintf('%d시 정각 갱신\n',sTh);
            elseif (j ~= 1 && mod(sTm,0.5) == 0) % 30분 일 때
                clear Sat_ar
                %         [Sat_ar] = GetSatPosGLO_new(EphGlo,gs,deltat);
                [Sat_ar] = GetSatPosGLO_my(EphGLO,gs,deltat);
                %         fprintf('%d시 30분 갱신\n',sTh);
            end
            
            %% 공통위성 선택
            GPSSats = intersect(GpsQMbs(:,2), GpsQMrv(:,2));        % GPS Satellites
            NoGPSSats = length(GPSSats);                            % Number of GPS's Sats, 1 Epoch
            GLOSats = intersect(GloQMbs(:,2), GloQMrv(:,2));        % GLO Satellites
            NoGLOSats = length(GLOSats);                            % Number of GLO's Sats, 1 Epoch
            
            %% gs, 공통위성 수, Base GPS, Rover GPS, Base GLO, Rover GLO
            visiSat(j,1) = gs; visiSat(j,2) = NoGPSSats;
            visiSat(j,3) = GpsSatsbs; visiSat(j,4) = GpsSatsrv;
            visiSat(j,5) = GpsSatsbs; visiSat(j,4) = GpsSatsrv;     % Number of visible Sats
            visiSat(j,7) = GpsSatsbs; visiSat(j,8) = GpsSatsrv;
            
            %% Zenith Hydrostatic Delay based on GPT
            ZHD = TropGPTh(vec_site, gw, gs);
            
            %% GPS 기준위성 RS와 다른위성 OS 선별/ SatsEl - c1(gs), c2(prn), c3(el)
            [GPSSatsEl, GPSindxRS] = PickRSel(gs, GPSSats, EphGPS, TruePosBs_PP);  % : GPS RS 기준위성
            GPSRS = GPSSats(GPSindxRS); RefGPSSV(j,1) = GPSRS;
            
            %% GLO 기준위성 RS와 다른위성 OS 선별/ SatsEl - c1(gs), c2(prn), c3(el)
            [GLOSatsEl, GLOindxRS] = PickRSelGLO(gs, LeapSec, GLOSats, EphGLO, Sat_ar, TauC, vec_site);  % : RS 기준위성
            GLORS = GLOSats(GLOindxRS); RefGLOSV(j,1) = GLORS;
            
            %% GPS 기준위성 좌표 계산 - Bs 기준
            icol = PickEPH(EphGPS, GPSRS ,gs);
            STT = GetSTTbrdc(gs, GPSRS, EphGPS, TruePosBs_PP);                 % 신호전달시간 계산
            tc = gs - STT;
            vec_GPSRS = GetSatPosNC(EphGPS, icol, tc);
            vec_GPSRS = RotSatPos(vec_GPSRS, STT);                                % 지구 자전효과 고려
            STT_rv = GetSTTbrdc(gs, GPSRS, EphGPS, x(1:3)');                 % 신호전달시간 계산
            tc_rv = gs - STT_rv;
            vec_GPSRS_rv = GetSatPosNC(EphGPS, icol, tc_rv);
            vec_GPSRS_rv = RotSatPos(vec_GPSRS_rv, STT_rv);                                % 지구 자전효과 고려
            
            %% GLO 기준위성 좌표 계산 - Bs 기준
            tc = gs - LeapSec;
            icol=PickEPH_GLO2(EphGLO, GLORS, tc);
            STT = GetSTTbrdcGLO2(Sat_ar,gs,GLORS, x(1:3));           % 신호전달시간 계산
            % LeapSecond & 신호전달 시간을 보정한 위성 위치 산출
            [vec_GLORS, vel_GLORS] = SatPosLS_STT(Sat_ar,gs,GLORS,LeapSec,STT,TauC);
            vec_GLORS = RotSatPos(vec_GLORS,STT);                                 % 지구 자전효과 고려
            
                SNRbsGPSRS = GpsQM_snrbs(find(GpsQM_snrbs(:,2) == GPSRS), 4);      % BASE GPS RS SNR 
                SNRrvGPSRS = GpsQM_snrrv(find(GpsQM_snrrv(:,2) == GPSRS), 4);      % ROVER GPS RS SNR
                SNRbsGLORS = GloQM_snrbs(find(GloQM_snrbs(:,2) == GLORS), 4);      % BASE GPS RS SNR 
                SNRrvGLORS = GloQM_snrrv(find(GloQM_snrrv(:,2) == GLORS), 4);      % ROVER GPS RS SNR
            %% Iterarion 시작
            for Iter = 1: MaxIter
                HTH = zeros(3,3);
                HTy = zeros(3,1);
                NoGPSSatsUsed = NoGPSSats;
                NoGLOSatsUsed = NoGLOSats;
                usedSatCount = 0;
                
                %% 각 위성에 대한 관측치, 계산치, H행렬 계산

                %% GPS DD part
                for kS = 1: NoGPSSats
                    if kS == GPSindxRS || GPSSatsEl(kS, 3) < eleCut
                        if GPSSatsEl(kS, 3) < eleCut
                            NoGPSSatsUsed = NoGPSSatsUsed - 1;
                            disp([GPSRS GPSSatsEl(kS, 2)])
                        end
                        continue
                    end
                    cntgps = cntgps + 1;
                    GPSOS = GPSSats(kS);                                          % Other GPS Sat PRN
                    OtherGPSSats(kS,1) = GPSOS;
                    SNRbsGPSOS = GpsQM_snrbs(find(GpsQM_snrbs(:,2) == GPSOS), 4);       % BASE GPS OS SNR matrix
                    SNRrvGPSOS = GpsQM_snrrv(find(GpsQM_snrrv(:,2) == GPSOS), 4);       % ROVER GPS OS SNR matrix
                    obs_BsRS = GpsQMbs(find(GpsQMbs(:, 2) == GPSRS), 4);
                    obs_RvRS = GpsQMrv(find(GpsQMrv(:, 2) == GPSRS), 4);
                    obs_BsOS = GpsQMbs(find(GpsQMbs(:, 2) == GPSOS), 4);
                    obs_RvOS = GpsQMrv(find(GpsQMrv(:, 2) == GPSOS), 4);
                    obs = (obs_BsRS - obs_RvRS) - (obs_BsOS - obs_RvOS);
                    
                    icol = PickEPH(EphGPS, GPSOS, gs);
                    STT = GetSTTbrdc(gs, GPSOS, EphGPS, x(1:3)'); % :  OS 위성위치를 추정좌표 기준으로 변경 11/9/14
%                     STT = GetSTTbrdc(gs, GPSOS, EphGPS, TruePosBs_PP);           % 신호전달시간 계산
                    tc = gs - STT;
                    vec_GPSOS = GetSatPosNC(EphGPS, icol, tc);
                    vec_GPSOS = RotSatPos(vec_GPSOS, STT);
                    STT_bs = GetSTTbrdc(gs, GPSOS, EphGPS, TruePosBs_PP);           % 신호전달시간 계산
                    tc_bs = gs - STT_bs;
                    vec_GPSOS_bs = GetSatPosNC(EphGPS, icol, tc_bs);
                    vec_GPSOS_bs = RotSatPos(vec_GPSOS_bs, STT_bs);
                    
                    %% DD 계산치 생성 파트 - 거리계산치를 각각 계산한 다음 DD 계산치 계산
                    vec_BsRS = vec_GPSRS - TruePosBs_PP';  com_BsRS = norm(vec_BsRS);
                    vec_RvRS = vec_GPSRS_rv - x(1:3);    com_RvRS = norm(vec_RvRS);
                    vec_BsOS = vec_GPSOS_bs - TruePosBs_PP';  com_BsOS = norm(vec_BsOS);
                    vec_RvOS = vec_GPSOS - x(1:3);    com_RvOS = norm(vec_RvOS);
                    com = (com_BsRS - com_RvRS) - (com_BsOS - com_RvOS);
                    y = obs -com;
                    
                    %% GPS RS, OS Elevation angle 저장
                    [azBsRS,elBsRS] = xyz2azel(vec_BsRS, AppLatBs, AppLonBs);
                    azelGPSRSbs(j,:) = [azBsRS, elBsRS];                        % GPS RS's Azimuth, Elevation Angle 저장
                    [azBsOS,elBsOS] = xyz2azel(vec_BsOS, AppLatBs, AppLonBs);
                    azelGPSOSbs(cntgps,:) = [azBsOS, elBsOS];                   % GPS OS's Azimuth, Elevation Angle 저장
                    
                    %% 가중치 파트
                    SNRGPSRS(cntgps,:) = [SNRbsGPSRS, SNRrvGPSRS,gs,GPSRS];     % GPS RS's Base, Rover SNR
                    SNRGPSOS(cntgps,:) = [SNRbsGPSOS, SNRrvGPSOS,gs,GPSOS];     % GPS OS's Base, Rover SNR
                    DDGPSel(cntgps,:) = [elBsRS, elBsOS,gs,GPSOS];
                    
%                     W =1;
                    W = DDMakeW_elsnr(SNRGPSRS(cntgps,:),SNRGPSOS(cntgps,:),DDGPSel(cntgps,:));
                    %% H 행렬 계산 파트
                    H(1,1) = vec_RvRS(1)/com_RvRS - vec_RvOS(1)/com_RvOS;
                    H(1,2) = vec_RvRS(2)/com_RvRS - vec_RvOS(2)/com_RvOS;
                    H(1,3) = vec_RvRS(3)/com_RvRS - vec_RvOS(3)/com_RvOS;
                    
                    HTH = HTH + H'*W*H;
                    HTy = HTy + H'*W*y;
                end
                
                
                %             OTHER{j,1} = OtherSats;
                % %             OTHER{j,2} = indxUsedSat;
                xhat = inv(HTH) * HTy;
                x = x + xhat;
                
                if norm(xhat) < EpsStop;
                    nEst = nEst + 1;
                    estm(nEst,1) =gs;
                    estm(nEst,2:4) =x;
                    estm(nEst,5) = NoGPSSats + NoGLOSats;
                    estm(nEst,6) = NoGPSSatsUsed + NoGLOSatsUsed;  % : 기준위성을 포함해야 함
                    estm(nEst,7) = 0;  % : snr issue 가 없을때
                    break;
                end
            end
            %% rover's longi, lati
            rover_gd(nEst,1) = gs;
            rover_gd(nEst,2:4) = xyz2gd(estm(nEst,2:4)); % rover의 xyz값 gd 로 변환
            AppLat = rover_gd(nEst,2); AppLon = rover_gd(nEst,3);
        end
        
    case 2
        
        for j = 1:NoEpochs
            % for j = 1:1
            gs = FinalTTs(j);
            idxQMbs = find(FinalQMbs(:,1) == gs);
            idxQMrv = find(FinalQMrv(:,1) == gs);
            QM_c1bs = FinalQMbs(idxQMbs,:);             % Base GPS/GLO C1, 1 Epoch
            QM_snrbs = FinalQM_snrbs(idxQMbs,:);        % Base GPS/GLO SNR, 1 Epoch
            QM_dopbs = FinalQM_dopbs(idxQMbs,:);        % Base GPS/GLO DOP, 1 Epoch
            QM_c1rv = FinalQMrv(idxQMrv,:);             % Rover GPS/GLO C1, 1 Epoch
            QM_snrrv = FinalQM_snrrv(idxQMrv,:);        % Rover GPS/GLO SNR, 1 Epoch
            QM_doprv = FinalQM_doprv(idxQMrv,:);        % Rover GPS/GLO DOP, 1 Epoch
            
            GpsQMbs = QM_c1bs(find(QM_c1bs(:,3) == 120),:); GpsSatsbs = length(GpsQMbs(:,2));               % Base GPS C1
            GloQMbs = QM_c1bs(find(QM_c1bs(:,3) == 320),:); GloSatsbs = length(GloQMbs(:,2));               % Base GLO C1
            GpsQM_snrbs = QM_snrbs(find(QM_snrbs(:,3) == 141),:); GpsSatsbs = length(GpsQM_snrbs(:,4));     % Base GPS SNR
            GloQM_snrbs = QM_snrbs(find(QM_snrbs(:,3) == 341),:); GloSatsbs = length(GloQM_snrbs(:,4));     % Base GLO SNR
            
            GpsQMrv = QM_c1rv(find(QM_c1rv(:,3) == 120),:); GpsSatsrv = length(GpsQMrv(:,2));               % Rover GPS C1
            GloQMrv = QM_c1rv(find(QM_c1rv(:,3) == 320),:); GloSatsrv = length(GloQMrv(:,2));               % Rover GLO C1
            GpsQM_snrrv = QM_snrrv(find(QM_snrrv(:,3) == 141),:); GpsSatsrv = length(GpsQM_snrrv(:,4));     % Rover GPS SNR
            GloQM_snrrv = QM_snrrv(find(QM_snrrv(:,3) == 341),:); GloSatsrv = length(GloQM_snrrv(:,4));     % Rover GLO SNR
            
            %% 해당 시각 Bs의 위치 값 찾기
            Base(j,:) = Bs(find(Bs(:,1) == gs),:);
            TruePosBs_PP = Base(j,2:4);
            base_gd(j,1) = gs;
            base_gd(j,2:4) = xyz2gd(TruePosBs_PP(1:3)); % Base의 xyz값 gd 로 변환
            AppLatBs = base_gd(j,2); AppLonBs = base_gd(j,3);
            vec_site = x(1:3)';
            
            %% GLO 위성 좌표 계산을 위한 적분 초기값 생성
            sT= mod(gs,86400)/3600;
            sTh = floor(sT); sTm = sT - sTh;
            if j == 1 % 시작시간 위성 위치 array
                %         [Sat_ar] = GetSatPosGLO_new(EphGlo,gs,deltat);
                [Sat_ar] = GetSatPosGLO_my(EphGLO,gs,deltat);
                %         fprintf('시작 시간 위치 계산\n');
            elseif (j ~= 1 && sTm == 0) % 정각일 때
                clear Sat_ar
                %         [Sat_ar] = GetSatPosGLO_new(EphGlo,gs,deltat);
                [Sat_ar] = GetSatPosGLO_my(EphGLO,gs,deltat);
                %         fprintf('%d시 정각 갱신\n',sTh);
            elseif (j ~= 1 && mod(sTm,0.5) == 0) % 30분 일 때
                clear Sat_ar
                %         [Sat_ar] = GetSatPosGLO_new(EphGlo,gs,deltat);
                [Sat_ar] = GetSatPosGLO_my(EphGLO,gs,deltat);
                %         fprintf('%d시 30분 갱신\n',sTh);
            end
            
            %% 공통위성 선택
            GPSSats = intersect(GpsQMbs(:,2), GpsQMrv(:,2));        % GPS Satellites
            NoGPSSats = length(GPSSats);                            % Number of GPS's Sats, 1 Epoch
            GLOSats = intersect(GloQMbs(:,2), GloQMrv(:,2));        % GLO Satellites
            NoGLOSats = length(GLOSats);                            % Number of GLO's Sats, 1 Epoch
            
            %% gs, 공통위성 수, Base GPS, Rover GPS, Base GLO, Rover GLO
            visiSat(j,1) = gs; visiSat(j,2) = NoGPSSats + NoGLOSats;
            visiSat(j,3) = GpsSatsbs; visiSat(j,4) = GpsSatsrv;
            visiSat(j,5) = GloSatsbs; visiSat(j,6) = GloSatsrv;     % Number of visible Sats
            visiSat(j,7) = GpsSatsbs +GloSatsbs; visiSat(j,8) = GpsSatsrv + GloSatsrv;
            
            %% Zenith Hydrostatic Delay based on GPT
            ZHD = TropGPTh(vec_site, gw, gs);
            
            %% GPS 기준위성 RS와 다른위성 OS 선별/ SatsEl - c1(gs), c2(prn), c3(el)
            [GPSSatsEl, GPSindxRS] = PickRSel(gs, GPSSats, EphGPS, TruePosBs_PP);  % : GPS RS 기준위성
            GPSRS = GPSSats(GPSindxRS); RefGPSSV(j,1) = GPSRS;
            
            %% GLO 기준위성 RS와 다른위성 OS 선별/ SatsEl - c1(gs), c2(prn), c3(el)
            [GLOSatsEl, GLOindxRS] = PickRSelGLO(gs, LeapSec, GLOSats, EphGLO, Sat_ar, TauC, vec_site);  % : RS 기준위성
            GLORS = GLOSats(GLOindxRS); RefGLOSV(j,1) = GLORS;
            
            %% GPS 기준위성 좌표 계산 - Bs 기준
            icol = PickEPH(EphGPS, GPSRS ,gs);
            STT = GetSTTbrdc(gs, GPSRS, EphGPS, TruePosBs_PP);                 % 신호전달시간 계산
            tc = gs - STT;
            vec_GPSRS = GetSatPosNC(EphGPS, icol, tc);
            vec_GPSRS = RotSatPos(vec_GPSRS, STT);                                % 지구 자전효과 고려
            
            %% GLO 기준위성 좌표 계산 - Bs 기준
            tc = gs - LeapSec;
            icol=PickEPH_GLO2(EphGLO, GLORS, tc);
            STT = GetSTTbrdcGLO2(Sat_ar,gs,GLORS, x(1:3));           % 신호전달시간 계산
            % LeapSecond & 신호전달 시간을 보정한 위성 위치 산출
            [vec_GLORS, vel_GLORS] = SatPosLS_STT(Sat_ar,gs,GLORS,LeapSec,STT,TauC);
            vec_GLORS = RotSatPos(vec_GLORS,STT);                                 % 지구 자전효과 고려
            
                SNRbsGPSRS = GpsQM_snrbs(find(GpsQM_snrbs(:,2) == GPSRS), 4);      % BASE GPS RS SNR 
                SNRrvGPSRS = GpsQM_snrrv(find(GpsQM_snrrv(:,2) == GPSRS), 4);      % ROVER GPS RS SNR
                SNRbsGLORS = GloQM_snrbs(find(GloQM_snrbs(:,2) == GLORS), 4);      % BASE GPS RS SNR 
                SNRrvGLORS = GloQM_snrrv(find(GloQM_snrrv(:,2) == GLORS), 4);      % ROVER GPS RS SNR
            %% Iterarion 시작
            for Iter = 1: MaxIter
                HTH = zeros(3,3);
                HTy = zeros(3,1);
                NoGPSSatsUsed = NoGPSSats;
                NoGLOSatsUsed = NoGLOSats;
                usedSatCount = 0;
                
                %% 각 위성에 대한 관측치, 계산치, H행렬 계산
                
                %% GLO DD part
                for kS = 1: NoGLOSats
                    if kS == GLOindxRS || GLOSatsEl(kS, 3) < eleCut
                        if GLOSatsEl(kS, 3) < eleCut
                            NoGLOSatsUsed = NoGLOSatsUsed - 1;
                            disp([GLORS GLOSatsEl(kS, 2)])
                        end
                        continue
                    end
                    cntglo = cntglo + 1;
                    GLOOS = GLOSats(kS);                                          % Other GLO Sat PRN
                    OtherGLOSats(kS,1) = GLOOS;
                    SNRbsGLOOS = GloQM_snrbs(find(GloQM_snrbs(:,2) == GLOOS), 4);       % BASE GLO OS SNR matrix
                    SNRrvGLOOS = GloQM_snrrv(find(GloQM_snrrv(:,2) == GLOOS), 4);       % ROVER GLO OS SNR matrix
                    
                    obs_BsRS = GloQMbs(find(GloQMbs(:, 2) == GLORS), 4);
                    obs_RvRS = GloQMrv(find(GloQMrv(:, 2) == GLORS), 4);
                    obs_BsOS = GloQMbs(find(GloQMbs(:, 2) == GLOOS), 4);
                    obs_RvOS = GloQMrv(find(GloQMrv(:, 2) == GLOOS), 4);
                    obs = (obs_BsRS - obs_RvRS) - (obs_BsOS - obs_RvOS);
                    tc = gs - LeapSec;
                    icol=PickEPH_GLO2(EphGLO, GLOOS, tc);
                    STT = GetSTTbrdcGLO2(Sat_ar,gs,GLOOS, x(1:3));           % 신호전달시간 계산

                    % LeapSecond & 신호전달 시간을 보정한 위성 위치 산출
                    [vec_GLOOS, vel_GLOOS] = SatPosLS_STT(Sat_ar,gs,GLOOS,LeapSec,STT,TauC);
                    vec_GLOOS = RotSatPos(vec_GLOOS,STT);                                 % 지구 자전효과 고려
                    
                    %% DD 계산치 생성 파트 - 거리계산치를 각각 계산한 다음 DD 계산치 계산
                    vec_BsRS = vec_GLORS' - TruePosBs_PP;  com_BsRS = norm(vec_BsRS);
                    vec_RvRS = vec_GLORS - x(1:3);    com_RvRS = norm(vec_RvRS);
                    vec_BsOS = vec_GLOOS' - TruePosBs_PP;  com_BsOS = norm(vec_BsOS);
                    vec_RvOS = vec_GLOOS - x(1:3);    com_RvOS = norm(vec_RvOS);
                    com = (com_BsRS - com_RvRS) - (com_BsOS - com_RvOS);
                    y = obs -com;
                    
                    %% GPS RS, OS Elevation angle 저장
                    [azBsRS,elBsRS] = xyz2azel(vec_BsRS, AppLatBs, AppLonBs);
                    azelGLORSbs(j,:) = [azBsRS, elBsRS];                        % GLO RS's Azimuth, Elevation Angle 저장
                    [azBsOS,elBsOS] = xyz2azel(vec_BsOS, AppLatBs, AppLonBs);
                    azelGLOOSbs(cntglo,:) = [azBsOS, elBsOS];                   % GLO OS's Azimuth, Elevation Angle 저장
                    
                    %% 가중치 파트
                    SNRGLORS(cntglo,:) = [SNRbsGLORS, SNRrvGLORS,gs,GLORS];     % GPS RS's Base, Rover SNR
                    SNRGLOOS(cntglo,:) = [SNRbsGLOOS, SNRrvGLOOS,gs,GLOOS];     % GPS OS's Base, Rover SNR
                    DDGLOel(cntglo,:) = [elBsRS, elBsOS,gs,GLOOS];
                    
                    W =1;
                    %% H 행렬 계산 파트
                    H(1,1) = vec_RvRS(1)/com_RvRS - vec_RvOS(1)/com_RvOS;
                    H(1,2) = vec_RvRS(2)/com_RvRS - vec_RvOS(2)/com_RvOS;
                    H(1,3) = vec_RvRS(3)/com_RvRS - vec_RvOS(3)/com_RvOS;
                    
                    HTH = HTH + H'*W*H;
                    HTy = HTy + H'*W*y;
                end
                %% GPS DD part
%                 for kS = 1: NoGPSSats
%                     if kS == GPSindxRS || GPSSatsEl(kS, 3) < eleCut
%                         if GPSSatsEl(kS, 3) < eleCut
%                             NoGPSSatsUsed = NoGPSSatsUsed - 1;
%                             disp([GPSRS GPSSatsEl(kS, 2)])
%                         end
%                         continue
%                     end
%                     cntgps = cntgps + 1;
%                     GPSOS = GPSSats(kS);                                          % Other GPS Sat PRN
%                     OtherGPSSats(kS,1) = GPSOS;
%                     SNRbsGPSOS = GpsQM_snrbs(find(GpsQM_snrbs(:,2) == GPSOS), 4);       % BASE GPS OS SNR matrix
%                     SNRrvGPSOS = GpsQM_snrrv(find(GpsQM_snrrv(:,2) == GPSOS), 4);       % ROVER GPS OS SNR matrix
%                     obs_BsRS = GpsQMbs(find(GpsQMbs(:, 2) == GPSRS), 4);
%                     obs_RvRS = GpsQMrv(find(GpsQMrv(:, 2) == GPSRS), 4);
%                     obs_BsOS = GpsQMbs(find(GpsQMbs(:, 2) == GPSOS), 4);
%                     obs_RvOS = GpsQMrv(find(GpsQMrv(:, 2) == GPSOS), 4);
%                     obs = (obs_BsRS - obs_RvRS) - (obs_BsOS - obs_RvOS);
%                     
%                     icol = PickEPH(EphGPS, GPSOS, gs);
%                     STT = GetSTTbrdc(gs, GPSOS, EphGPS, x(1:3)'); % :  OS 위성위치를 추정좌표 기준으로 변경 11/9/14
%                     tc = gs - STT;
%                     vec_GPSOS = GetSatPosNC(EphGPS, icol, tc);
%                     vec_GPSOS = RotSatPos(vec_GPSOS, STT);
%                     
%                     %% DD 계산치 생성 파트 - 거리계산치를 각각 계산한 다음 DD 계산치 계산
%                     vec_BsRS = vec_GPSRS - TruePosBs_PP;  com_BsRS = norm(vec_BsRS);
%                     vec_RvRS = vec_GPSRS - x(1:3)';    com_RvRS = norm(vec_RvRS);
%                     vec_BsOS = vec_GPSOS - TruePosBs_PP;  com_BsOS = norm(vec_BsOS);
%                     vec_RvOS = vec_GPSOS - x(1:3)';    com_RvOS = norm(vec_RvOS);
%                     com = (com_BsRS - com_RvRS) - (com_BsOS - com_RvOS);
%                     y = obs -com;
%                     
%                     %% GPS RS, OS Elevation angle 저장
%                     [azBsRS,elBsRS] = xyz2azel(vec_BsRS, AppLatBs, AppLonBs);
%                     azelGPSRSbs(j,:) = [azBsRS, elBsRS];                        % GPS RS's Azimuth, Elevation Angle 저장
%                     [azBsOS,elBsOS] = xyz2azel(vec_BsOS, AppLatBs, AppLonBs);
%                     azelGPSOSbs(cntgps,:) = [azBsOS, elBsOS];                   % GPS OS's Azimuth, Elevation Angle 저장
%                     
%                     %% 가중치 파트
%                     SNRGPSRS(cntgps,:) = [SNRbsGPSRS, SNRrvGPSRS,gs,GPSRS];     % GPS RS's Base, Rover SNR
%                     SNRGPSOS(cntgps,:) = [SNRbsGPSOS, SNRrvGPSOS,gs,GPSOS];     % GPS OS's Base, Rover SNR
%                     DDGPSel(cntgps,:) = [elBsRS, elBsOS,gs,GPSOS];
%                     
%                     W =1;
%                     %% H 행렬 계산 파트
%                     H(1,1) = vec_RvRS(1)/com_RvRS - vec_RvOS(1)/com_RvOS;
%                     H(1,2) = vec_RvRS(2)/com_RvRS - vec_RvOS(2)/com_RvOS;
%                     H(1,3) = vec_RvRS(3)/com_RvRS - vec_RvOS(3)/com_RvOS;
%                     
%                     HTH = HTH + H'*W*H;
%                     HTy = HTy + H'*W*y;
%                 end

                xhat = inv(HTH) * HTy;
                x = x + xhat;
                
                if norm(xhat) < EpsStop;
                    nEst = nEst + 1;
                    estm(nEst,1) =gs;
                    estm(nEst,2:4) =x;
                    estm(nEst,5) = NoGPSSats + NoGLOSats;
                    estm(nEst,6) = NoGPSSatsUsed + NoGLOSatsUsed;  % : 기준위성을 포함해야 함
                    estm(nEst,7) = 0;  % : snr issue 가 없을때
                    break;
                end
            end
            %% rover's longi, lati
            rover_gd(nEst,1) = gs;
            rover_gd(nEst,2:4) = xyz2gd(estm(nEst,2:4)); % rover의 xyz값 gd 로 변환
            AppLat = rover_gd(nEst,2); AppLon = rover_gd(nEst,3);
        end
        
    case 3
        
        for j = 1:NoEpochs
            % for j = 1:1
            gs = FinalTTs(j);
            idxQMbs = find(FinalQMbs(:,1) == gs);
            idxQMrv = find(FinalQMrv(:,1) == gs);
            QM_c1bs = FinalQMbs(idxQMbs,:);             % Base GPS/GLO C1, 1 Epoch
            QM_snrbs = FinalQM_snrbs(idxQMbs,:);        % Base GPS/GLO SNR, 1 Epoch
            QM_dopbs = FinalQM_dopbs(idxQMbs,:);        % Base GPS/GLO DOP, 1 Epoch
            QM_c1rv = FinalQMrv(idxQMrv,:);             % Rover GPS/GLO C1, 1 Epoch
            QM_snrrv = FinalQM_snrrv(idxQMrv,:);        % Rover GPS/GLO SNR, 1 Epoch
            QM_doprv = FinalQM_doprv(idxQMrv,:);        % Rover GPS/GLO DOP, 1 Epoch
            
            GpsQMbs = QM_c1bs(find(QM_c1bs(:,3) == 120),:); GpsSatsbs = length(GpsQMbs(:,2));               % Base GPS C1
            GloQMbs = QM_c1bs(find(QM_c1bs(:,3) == 320),:); GloSatsbs = length(GloQMbs(:,2));               % Base GLO C1
            GpsQM_snrbs = QM_snrbs(find(QM_snrbs(:,3) == 141),:); GpsSatsbs = length(GpsQM_snrbs(:,4));     % Base GPS SNR
            GloQM_snrbs = QM_snrbs(find(QM_snrbs(:,3) == 341),:); GloSatsbs = length(GloQM_snrbs(:,4));     % Base GLO SNR
            
            GpsQMrv = QM_c1rv(find(QM_c1rv(:,3) == 120),:); GpsSatsrv = length(GpsQMrv(:,2));               % Rover GPS C1
            GloQMrv = QM_c1rv(find(QM_c1rv(:,3) == 320),:); GloSatsrv = length(GloQMrv(:,2));               % Rover GLO C1
            GpsQM_snrrv = QM_snrrv(find(QM_snrrv(:,3) == 141),:); GpsSatsrv = length(GpsQM_snrrv(:,4));     % Rover GPS SNR
            GloQM_snrrv = QM_snrrv(find(QM_snrrv(:,3) == 341),:); GloSatsrv = length(GloQM_snrrv(:,4));     % Rover GLO SNR
            
            %% 해당 시각 Bs의 위치 값 찾기
            Base(j,:) = Bs(find(Bs(:,1) == gs),:);
            TruePosBs_PP = Base(j,2:4);
            base_gd(j,1) = gs;
            base_gd(j,2:4) = xyz2gd(TruePosBs_PP(1:3)); % Base의 xyz값 gd 로 변환
            AppLatBs = base_gd(j,2); AppLonBs = base_gd(j,3);
            vec_site = x(1:3)';
            
            %% GLO 위성 좌표 계산을 위한 적분 초기값 생성
            sT= mod(gs,86400)/3600;
            sTh = floor(sT); sTm = sT - sTh;
            if j == 1 % 시작시간 위성 위치 array
                %         [Sat_ar] = GetSatPosGLO_new(EphGlo,gs,deltat);
                [Sat_ar] = GetSatPosGLO_my(EphGLO,gs,deltat);
                %         fprintf('시작 시간 위치 계산\n');
            elseif (j ~= 1 && sTm == 0) % 정각일 때
                clear Sat_ar
                %         [Sat_ar] = GetSatPosGLO_new(EphGlo,gs,deltat);
                [Sat_ar] = GetSatPosGLO_my(EphGLO,gs,deltat);
                %         fprintf('%d시 정각 갱신\n',sTh);
            elseif (j ~= 1 && mod(sTm,0.5) == 0) % 30분 일 때
                clear Sat_ar
                %         [Sat_ar] = GetSatPosGLO_new(EphGlo,gs,deltat);
                [Sat_ar] = GetSatPosGLO_my(EphGLO,gs,deltat);
                %         fprintf('%d시 30분 갱신\n',sTh);
            end
            
            %% 공통위성 선택
            GPSSats = intersect(GpsQMbs(:,2), GpsQMrv(:,2));        % GPS Satellites
            NoGPSSats = length(GPSSats);                            % Number of GPS's Sats, 1 Epoch
            GLOSats = intersect(GloQMbs(:,2), GloQMrv(:,2));        % GLO Satellites
            NoGLOSats = length(GLOSats);                            % Number of GLO's Sats, 1 Epoch
            
            %% gs, 공통위성 수, Base GPS, Rover GPS, Base GLO, Rover GLO
            visiSat(j,1) = gs; visiSat(j,2) = NoGPSSats + NoGLOSats;
            visiSat(j,3) = GpsSatsbs; visiSat(j,4) = GpsSatsrv;
            visiSat(j,5) = GloSatsbs; visiSat(j,6) = GloSatsrv;     % Number of visible Sats
            visiSat(j,7) = GpsSatsbs +GloSatsbs; visiSat(j,8) = GpsSatsrv + GloSatsrv;
            
            %% Zenith Hydrostatic Delay based on GPT
            ZHD = TropGPTh(vec_site, gw, gs);
            
            %% GPS 기준위성 RS와 다른위성 OS 선별/ SatsEl - c1(gs), c2(prn), c3(el)
            [GPSSatsEl, GPSindxRS] = PickRSel(gs, GPSSats, EphGPS, TruePosBs_PP);  % : GPS RS 기준위성
            GPSRS = GPSSats(GPSindxRS); RefGPSSV(j,1) = GPSRS;
            
            %% GLO 기준위성 RS와 다른위성 OS 선별/ SatsEl - c1(gs), c2(prn), c3(el)
            [GLOSatsEl, GLOindxRS] = PickRSelGLO(gs, LeapSec, GLOSats, EphGLO, Sat_ar, TauC, vec_site);  % : RS 기준위성
            GLORS = GLOSats(GLOindxRS); RefGLOSV(j,1) = GLORS;
            
            %% GPS 기준위성 좌표 계산 - Bs 기준
            icol = PickEPH(EphGPS, GPSRS ,gs);
            STT = GetSTTbrdc(gs, GPSRS, EphGPS, TruePosBs_PP);                 % 신호전달시간 계산
            tc = gs - STT;
            vec_GPSRS = GetSatPosNC(EphGPS, icol, tc);
            vec_GPSRS = RotSatPos(vec_GPSRS, STT);                                % 지구 자전효과 고려
            
            %% GLO 기준위성 좌표 계산 - Bs 기준
            tc = gs - LeapSec;
            icol=PickEPH_GLO2(EphGLO, GLORS, tc);
            STT = GetSTTbrdcGLO2(Sat_ar,gs,GLORS, x(1:3));           % 신호전달시간 계산
            % LeapSecond & 신호전달 시간을 보정한 위성 위치 산출
            [vec_GLORS, vel_GLORS] = SatPosLS_STT(Sat_ar,gs,GLORS,LeapSec,STT,TauC);
            vec_GLORS = RotSatPos(vec_GLORS,STT);                                 % 지구 자전효과 고려
            
                SNRbsGPSRS = GpsQM_snrbs(find(GpsQM_snrbs(:,2) == GPSRS), 4);      % BASE GPS RS SNR 
                SNRrvGPSRS = GpsQM_snrrv(find(GpsQM_snrrv(:,2) == GPSRS), 4);      % ROVER GPS RS SNR
                SNRbsGLORS = GloQM_snrbs(find(GloQM_snrbs(:,2) == GLORS), 4);      % BASE GPS RS SNR 
                SNRrvGLORS = GloQM_snrrv(find(GloQM_snrrv(:,2) == GLORS), 4);      % ROVER GPS RS SNR
            %% Iterarion 시작
            for Iter = 1: MaxIter
                HTH = zeros(3,3);
                HTy = zeros(3,1);
                NoGPSSatsUsed = NoGPSSats;
                NoGLOSatsUsed = NoGLOSats;
                usedSatCount = 0;
                
                %% 각 위성에 대한 관측치, 계산치, H행렬 계산
                
                %% GLO DD part
                for kS = 1: NoGLOSats
                    if kS == GLOindxRS || GLOSatsEl(kS, 3) < eleCut
                        if GLOSatsEl(kS, 3) < eleCut
                            NoGLOSatsUsed = NoGLOSatsUsed - 1;
                            disp([GLORS GLOSatsEl(kS, 2)])
                        end
                        continue
                    end
                    cntglo = cntglo + 1;
                    GLOOS = GLOSats(kS);                                          % Other GLO Sat PRN
                    OtherGLOSats(kS,1) = GLOOS;
                    SNRbsGLOOS = GloQM_snrbs(find(GloQM_snrbs(:,2) == GLOOS), 4);       % BASE GLO OS SNR matrix
                    SNRrvGLOOS = GloQM_snrrv(find(GloQM_snrrv(:,2) == GLOOS), 4);       % ROVER GLO OS SNR matrix
                    
                    obs_BsRS = GloQMbs(find(GloQMbs(:, 2) == GLORS), 4);
                    obs_RvRS = GloQMrv(find(GloQMrv(:, 2) == GLORS), 4);
                    obs_BsOS = GloQMbs(find(GloQMbs(:, 2) == GLOOS), 4);
                    obs_RvOS = GloQMrv(find(GloQMrv(:, 2) == GLOOS), 4);
                    obs = (obs_BsRS - obs_RvRS) - (obs_BsOS - obs_RvOS);
                    tc = gs - LeapSec;
                    icol=PickEPH_GLO2(EphGLO, GLOOS, tc);
                    STT = GetSTTbrdcGLO2(Sat_ar,gs,GLOOS, x(1:3));           % 신호전달시간 계산
                    % LeapSecond & 신호전달 시간을 보정한 위성 위치 산출
                    [vec_GLOOS, vel_GLOOS] = SatPosLS_STT(Sat_ar,gs,GLOOS,LeapSec,STT,TauC);
                    vec_GLOOS = RotSatPos(vec_GLOOS,STT);                                 % 지구 자전효과 고려
                    
                    %% DD 계산치 생성 파트 - 거리계산치를 각각 계산한 다음 DD 계산치 계산
                    vec_BsRS = vec_GLORS - TruePosBs_PP;  com_BsRS = norm(vec_BsRS);
                    vec_RvRS = vec_GLORS - x(1:3)';    com_RvRS = norm(vec_RvRS);
                    vec_BsOS = vec_GLOOS - TruePosBs_PP;  com_BsOS = norm(vec_BsOS);
                    vec_RvOS = vec_GLOOS - x(1:3)';    com_RvOS = norm(vec_RvOS);
                    com = (com_BsRS - com_RvRS) - (com_BsOS - com_RvOS);
                    y = obs -com;
                    
                    %% GPS RS, OS Elevation angle 저장
                    [azBsRS,elBsRS] = xyz2azel(vec_BsRS, AppLatBs, AppLonBs);
                    azelGLORSbs(j,:) = [azBsRS, elBsRS];                        % GLO RS's Azimuth, Elevation Angle 저장
                    [azBsOS,elBsOS] = xyz2azel(vec_BsOS, AppLatBs, AppLonBs);
                    azelGLOOSbs(cntglo,:) = [azBsOS, elBsOS];                   % GLO OS's Azimuth, Elevation Angle 저장
                    
                    %% 가중치 파트
                    SNRGLORS(cntglo,:) = [SNRbsGLORS, SNRrvGLORS,gs,GLORS];     % GPS RS's Base, Rover SNR
                    SNRGLOOS(cntglo,:) = [SNRbsGLOOS, SNRrvGLOOS,gs,GLOOS];     % GPS OS's Base, Rover SNR
                    DDGLOel(cntglo,:) = [elBsRS, elBsOS,gs,GLOOS];
                    
                    W =1;
                    %% H 행렬 계산 파트
                    H(1,1) = vec_RvRS(1)/com_RvRS - vec_RvOS(1)/com_RvOS;
                    H(1,2) = vec_RvRS(2)/com_RvRS - vec_RvOS(2)/com_RvOS;
                    H(1,3) = vec_RvRS(3)/com_RvRS - vec_RvOS(3)/com_RvOS;
                    
                    HTH = HTH + H'*W*H;
                    HTy = HTy + H'*W*y;
                end
                %% GPS DD part
                for kS = 1: NoGPSSats
                    if kS == GPSindxRS || GPSSatsEl(kS, 3) < eleCut
                        if GPSSatsEl(kS, 3) < eleCut
                            NoGPSSatsUsed = NoGPSSatsUsed - 1;
                            disp([GPSRS GPSSatsEl(kS, 2)])
                        end
                        continue
                    end
                    cntgps = cntgps + 1;
                    GPSOS = GPSSats(kS);                                          % Other GPS Sat PRN
                    OtherGPSSats(kS,1) = GPSOS;
                    SNRbsGPSOS = GpsQM_snrbs(find(GpsQM_snrbs(:,2) == GPSOS), 4);       % BASE GPS OS SNR matrix
                    SNRrvGPSOS = GpsQM_snrrv(find(GpsQM_snrrv(:,2) == GPSOS), 4);       % ROVER GPS OS SNR matrix
                    obs_BsRS = GpsQMbs(find(GpsQMbs(:, 2) == GPSRS), 4);
                    obs_RvRS = GpsQMrv(find(GpsQMrv(:, 2) == GPSRS), 4);
                    obs_BsOS = GpsQMbs(find(GpsQMbs(:, 2) == GPSOS), 4);
                    obs_RvOS = GpsQMrv(find(GpsQMrv(:, 2) == GPSOS), 4);
                    obs = (obs_BsRS - obs_RvRS) - (obs_BsOS - obs_RvOS);
                    
                    icol = PickEPH(EphGPS, GPSOS, gs);
                    STT = GetSTTbrdc(gs, GPSOS, EphGPS, x(1:3)'); % :  OS 위성위치를 추정좌표 기준으로 변경 11/9/14
                    tc = gs - STT;
                    vec_GPSOS = GetSatPosNC(EphGPS, icol, tc);
                    vec_GPSOS = RotSatPos(vec_GPSOS, STT);
                    
                    %% DD 계산치 생성 파트 - 거리계산치를 각각 계산한 다음 DD 계산치 계산
                    vec_BsRS = vec_GPSRS - TruePosBs_PP;  com_BsRS = norm(vec_BsRS);
                    vec_RvRS = vec_GPSRS - x(1:3)';    com_RvRS = norm(vec_RvRS);
                    vec_BsOS = vec_GPSOS - TruePosBs_PP;  com_BsOS = norm(vec_BsOS);
                    vec_RvOS = vec_GPSOS - x(1:3)';    com_RvOS = norm(vec_RvOS);
                    com = (com_BsRS - com_RvRS) - (com_BsOS - com_RvOS);
                    y = obs -com;
                    
                    %% GPS RS, OS Elevation angle 저장
                    [azBsRS,elBsRS] = xyz2azel(vec_BsRS, AppLatBs, AppLonBs);
                    azelGPSRSbs(j,:) = [azBsRS, elBsRS];                        % GPS RS's Azimuth, Elevation Angle 저장
                    [azBsOS,elBsOS] = xyz2azel(vec_BsOS, AppLatBs, AppLonBs);
                    azelGPSOSbs(cntgps,:) = [azBsOS, elBsOS];                   % GPS OS's Azimuth, Elevation Angle 저장
                    
                    %% 가중치 파트
                    SNRGPSRS(cntgps,:) = [SNRbsGPSRS, SNRrvGPSRS,gs,GPSRS];     % GPS RS's Base, Rover SNR
                    SNRGPSOS(cntgps,:) = [SNRbsGPSOS, SNRrvGPSOS,gs,GPSOS];     % GPS OS's Base, Rover SNR
                    DDGPSel(cntgps,:) = [elBsRS, elBsOS,gs,GPSOS];
                    
                    W =1;
                    %% H 행렬 계산 파트
                    H(1,1) = vec_RvRS(1)/com_RvRS - vec_RvOS(1)/com_RvOS;
                    H(1,2) = vec_RvRS(2)/com_RvRS - vec_RvOS(2)/com_RvOS;
                    H(1,3) = vec_RvRS(3)/com_RvRS - vec_RvOS(3)/com_RvOS;
                    
                    HTH = HTH + H'*W*H;
                    HTy = HTy + H'*W*y;
                end
                
                xhat = inv(HTH) * HTy;
                x = x + xhat;
                
                if norm(xhat) < EpsStop;
                    nEst = nEst + 1;
                    estm(nEst,1) =gs;
                    estm(nEst,2:4) =x;
                    estm(nEst,5) = NoGPSSats + NoGLOSats;
                    estm(nEst,6) = NoGPSSatsUsed + NoGLOSatsUsed;  % : 기준위성을 포함해야 함
                    estm(nEst,7) = 0;  % : snr issue 가 없을때
                    break;
                end
            end
            %% rover's longi, lati
            rover_gd(nEst,1) = gs;
            rover_gd(nEst,2:4) = xyz2gd(estm(nEst,2:4)); % rover의 xyz값 gd 로 변환
            AppLat = rover_gd(nEst,2); AppLon = rover_gd(nEst,3);
        end
end
% end
[DDdXYZ, DDdNEV, DDdis, DDrms, DDstd] = PostErrorsDDmv2(estm, Base, Truedis,0, 5, visiSat);
savefileloc = 'D:\ppsoln\ppsoln\matlab\data\DD\DDfigures\';
savefigurename = strcat(savefileloc,obsfileBs(1:8),'_',num2str(System));
savematname = strcat(savefileloc,obsfileBs(1:8),'_',num2str(System),'.mat');
saveas(gcf,savefigurename,'fig');
% saveas(gcf,savefigurename,'png');
save(savematname,'estm','Base','Truedis','visiSat');

toc


% 동일 epoch base, rover 구글 Plot
% figure(200)
% grid on
% hold on
% axis([min(base_gd(:,3))-0.0001 max(base_gd(:,3))+0.0001 min(base_gd(:,2))-0.0001 max(base_gd(:,2))+0.0001])
% plot_google_map;
% axis([min(base_gd(:,3))-0.0001 max(base_gd(:,3))+0.0001 min(base_gd(:,2))-0.0001 max(base_gd(:,2))+0.0001])
% finalepoch = intersect(base_gd(:,1), rover_gd(:,1));
% i = 1;
% for i = 1: length(finalepoch)
%     epoch = finalepoch(i);
%     Bs_gd = base_gd(find(base_gd(:,1) == epoch), 2:4);
%     Rv_gd = rover_gd(find(rover_gd(:,1) == epoch), 2:4);
%     setlon(i,:) = [Bs_gd(2), Rv_gd(2)];
%     setla(i,:) = [Bs_gd(1), Rv_gd(1)];
%     figure(200)
%     plot(setlon(i,:), setla(i,:),'r-');
%     grid on; hold on;
%     plot(Bs_gd(2), Bs_gd(1),'r.','MarkerSize',20)
%     plot(Rv_gd(2), Rv_gd(1),'b.','MarkerSize',20)
% %     legend('epoch','Left','Right')
% end