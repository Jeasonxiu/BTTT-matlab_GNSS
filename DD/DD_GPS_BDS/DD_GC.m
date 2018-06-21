%% 코드의사거리 이중차분 알고리즘
% 07/01/2016 : Joonseong
close all; clear all;
% load('DD_gc_test.mat');

%% 불변 변수 설정 : 빛의 속도, 관측치
CCC = 299792458.;   % CCC = Speed of Light [m/s]
Freq_L1 = 1575.42e6;        % L1 주파수
Lambda_L1 = CCC/Freq_L1;    % L1 파장
ObsType = 120;      % 사용할 관측지 설정 - 120 : C/A = C1
ObsType2 = 141;      % 사용할 관측치 설정 - 141: snr = S1
%% 임계고도각 설정
eleCut = 15;
%% 윤초 핸들링
LeapSecBDS = 14;
%% True Distance
% Truedis = 1.41;         % 신항대로 실험
Truedis = 9.9209;           % 대성 디폴리스 A,B

%% 계산 날짜 설정
DOY = 025; YY  = 17;
[gw, gd] = ydoy2gwgd(YY, DOY); %: GPS WEEK 결정

%% 항법메시지 호출
navfile = strcat('brdm', num2str(DOY,'%03d'), '0.', num2str(YY,'%02d'), 'p');   %: Navigation RINEX file
% [eph, trashPRN, trashT]=ReadEPH_all(navfile);

%% 사용할 관측치 선정
g_ObsType = 111; % gps C1
g_ObsType_snr = 141;
c_ObsType = 211; % bds C1
c_ObsType_snr = 241;

%% QM 파일 핸들링
QMfileBs = 'QM170125_A';
QMfileRv = 'QM170125_B';

%% QM 파일 읽어들여서 행렬로 저장하고, 사용할 관측치 추출
% 베이스 QM
[arrQM_Bs, FinalPRNs_Bs, FinalTTs_Bs] = ReadQM(QMfileBs);
arrQM_Bs = arrQM_Bs(find(arrQM_Bs(:,3) < 300),:);
QM_Bs = SelectQM_gc(arrQM_Bs, g_ObsType, c_ObsType);
QM_Bs_snr = SelectQM_gc(arrQM_Bs, g_ObsType_snr, c_ObsType_snr);
% 로버 QM
[arrQM_Rv, FinalPRNs_Rv, FinalTTs_Rv] = ReadQM(QMfileRv);
arrQM_Rv = arrQM_Rv(find(arrQM_Rv(:,3) < 300),:);
QM_Rv = SelectQM_gc(arrQM_Rv, g_ObsType, c_ObsType);
QM_Rv_snr = SelectQM_gc(arrQM_Rv, g_ObsType_snr, c_ObsType_snr);

%% 좌표 참값
% TruePos = [-3058799.61420451,4083265.35912516,3814946.87192938];
TruePos = [-3041235.578 4053941.677 3859881.013]; 

%% 항법메시지를 읽어들여서 행렬로 저장하고, Klobuchar 모델 추출
gps_nav = strcat(navfile(1:3),'c',navfile(5:8),'.',navfile(10:11),'n');
fid = fopen(gps_nav,'r');
if fid == -1
    al = zeros(4,1); be = zeros(4,1);
else
    [al, be] = GetALBE(navfile);
end

%% 라이넥스 파일에서 대략적인 관측소 좌표를 뽑아냄
AppPos = TruePos;

%% 라이넥스 파일에서 Base Station의 좌표를 계산함
% Bs = PP_gc(QMfileBs,eph,TruePos,DOY,YY);              % without Correction
load('DD_gc_Bs.mat');
%% 두 QM 파일에서 공통시각(epoch) 추출
FinalTTs = intersect(Bs(:, 1), QM_Rv(:, 1));


%% 추정에 필요한 초기치 설정
MaxIter = 10;
EpsStop = 1e-4;
x = [AppPos'; 1]
x_ = AppPos';


%% 추정과정 시작
NoEpochs = length(FinalTTs);
estm = zeros(NoEpochs, 4);  % : c1(gs), c2/3/4(X/Y/Z of Rv), c5(#Sats_obs), c6(#Sats_used)
nEst = 0;
cnt = 1;

for j = 1:NoEpochs
%     for j = 1:100
% for j = 100:500
    gs = FinalTTs(j);
    %% 해당 시각 Bs의 위치 값 찾기(PP-LS)
    Base(j,:) = Bs(find(Bs(:,1) == gs),:);
    TruePosBs_PP = Base(j,2:4);                                             % 베이스 좌표를 기준좌표로 설정
    %% base의 gd
    base_gd(j,1) = gs;
    base_gd(j,2:4) = xyz2gd(TruePosBs_PP(1:3));                             % Base의 xyz값 gd 로 변환
    AppLatBs = base_gd(j,2); AppLonBs = base_gd(j,3);
    %% 해당 시각 gs의 베이스 관측치 선별
    indexQM1 = find(QM_Bs(:,1) == gs);
    QM_Bs_1e = QM_Bs(indexQM1,:);                                           % Base Pseudo-Range(gs)
    QM_Bs_snr_1e = QM_Bs_snr(indexQM1,:);                                   % Base SNR(gs)
    
    %% 해당 시각 gs의 로버 관측치 선별
    indexQM2 = find(QM_Rv(:,1) == gs);
    QM_Rv_1e = QM_Rv(indexQM2,:);                                           % Rover Pseudo-Range(gs)
    QM_Rv_snr_1e = QM_Rv_snr(indexQM2,:);                                   % Rover SNR(gs)
         
    Sats = intersect(QM_Bs_1e(:, 2), QM_Rv_1e(:, 2));
    NoSats = length(Sats);
    
    %% 기준위성 RS와 다른위성 OS 선별/ SatsEl - c1(gs), c2(prn), c3(el)
    [SatsEl, indxRS] = PickRSel_gc(gs, Sats, eph, TruePosBs_PP);  % : RS 기준위성
    GPS = SatsEl(find(SatsEl(:,2) < 200),:);                                % SatsEl 중 GPS 위성군 선택
    BDS = SatsEl(find(SatsEl(:,2) > 200),:);                                % SatsEl 중 BDS 위성군 선택
    GPSRS = GPS(find(GPS(:,3) == max(GPS(:,3))),2);                         % GPS Reference Sat 선택
    BDSRS = BDS(find(BDS(:,3) == max(BDS(:,3))),2);                         % BDS Reference Sat 선택
    GPSindexRS = find(GPS(:,3) == max(GPS(:,3)));                           % GPS Reference Sat Index
    BDSindexRS = find(BDS(:,3) == max(BDS(:,3)));                           % BDS Reference Sat Index
    RefSV(j,1:2) = [GPSRS, BDSRS];                                          % 기준위성 PRN 저장
    
    %% gs, 공통위성 수, Base GPS, Rover GPS, Base GLO, Rover GLO...
    %% [gs, 전체 관측 위성수, 베이스 GPS 관측 수, 로버 GPS 관측수,...
    %%  베이스 BDS 관측 수, 로버 BDS 관측 수, 베이스 전체 관측수, 로버 전체 관측수]
    visiSat(j,1) = gs; visiSat(j,2) = length(GPS(:,1)) + length(BDS(:,1));
    visiSat(j,3) = length(QM_Bs_1e(find(QM_Bs_1e(:,2) < 200), 2)); 
    visiSat(j,4) = length(QM_Rv_1e(find(QM_Rv_1e(:,2) < 200), 2)); 
    visiSat(j,5) = length(QM_Bs_1e(find(QM_Bs_1e(:,2) > 200), 2)); 
    visiSat(j,6) = length(QM_Rv_1e(find(QM_Rv_1e(:,2) > 200), 2)); 
    visiSat(j,7) = visiSat(j,3)+visiSat(j,5); 
    visiSat(j,8) = visiSat(j,4)+visiSat(j,6); 
    
    %% GPS/BDS 기준 위성 snr 값 추출
    GPS_snr_Bs_RS = QM_Bs_snr_1e(find(QM_Bs_snr_1e(:,2) == GPSRS),4);           % 베이스, GPS 기준 위성 snr 값
    BDS_snr_Bs_RS = QM_Bs_snr_1e(find(QM_Bs_snr_1e(:,2) == BDSRS),4);           % 베이스, BDS 기준 위성 snr 값
    GPS_snr_Rv_RS = QM_Rv_snr_1e(find(QM_Rv_snr_1e(:,2) == GPSRS),4);           % 베이스, GPS 기준 위성 snr 값
    BDS_snr_Rv_RS = QM_Rv_snr_1e(find(QM_Rv_snr_1e(:,2) == BDSRS),4);           % 베이스, BDS 기준 위성 snr 값
    
    %% Iteration 시작
    for Iter = 1:MaxIter
        
%         HTH = zeros(3,3);
%         HTy = zeros(3,1);

        cnt2=1;
        %% GPS 기준위성 좌표 먼저 계산 - Bs 기준
        icol = PickEPH(eph, GPSRS ,gs);
        STT = GetSTTbrdm(gs, eph, icol, TruePosBs_PP');                     % 신호전달시간 계산
        tc = gs - STT;
        vec_GPSRS = GetSatPosNC_GC(eph, icol, tc);
        vec_GPSRS = RotSatPos(vec_GPSRS, STT);                              % 지구 자전효과 고려
        %% BDS 기준위성 좌표 먼저 계산 - Bs 기준
        icol = PickEPH(eph, BDSRS ,gs);
        STT = GetSTTbrdm(gs, eph, icol, TruePosBs_PP');                     % 신호전달시간 계산
        tc = gs - STT;
        vec_BDSRS = GetSatPosNC_GC(eph, icol, tc);
        vec_BDSRS = RotSatPos(vec_BDSRS, STT);                              % 지구 자전효과 고려
        NoGPSSatsUsed = length(GPS(:,1));
        NoBDSSatsUsed = length(BDS(:,1));
        
        %% GPS DD part
        for kS = 1:length(GPS(:,1))
            if kS == GPSindexRS || GPS(kS, 3) < eleCut
                if GPS(kS, 3) < eleCut
                    NoGPSSatsUsed = NoGPSSatsUsed - 1;
                    disp([GPSRS GPS(kS, 2)])
                end
                continue
            end
            %% DD 관측치 생성 파트
            GPSOS = GPS(kS,2);
            OtherGPSSats(kS,1) = GPSOS;
            GPS_snr_Bs_OS = QM_Bs_snr_1e(find(QM_Bs_snr_1e(:,2) == GPSOS),4);           % 베이스, GPS 기준 위성 snr 값
            GPS_snr_Rv_OS = QM_Rv_snr_1e(find(QM_Rv_snr_1e(:,2) == GPSOS),4);           % 베이스, GPS 기준 위성 snr 값
                        
            S1BsOS = QM_Bs_snr_1e(find(QM_Bs_snr_1e(:,2) == GPSOS), 4);      % BASE OS SNR matrix
            S1RvOS = QM_Rv_snr_1e(find(QM_Rv_snr_1e(:,2) == GPSOS), 4);         % ROVER OS SNR matrix
            
            obs_BsRS = QM_Bs_1e(find(QM_Bs_1e(:, 2) == GPSRS), 4);
            obs_RvRS = QM_Rv_1e(find(QM_Rv_1e(:, 2) == GPSRS), 4);
            obs_BsOS = QM_Bs_1e(find(QM_Bs_1e(:, 2) == GPSOS), 4);
            obs_RvOS = QM_Rv_1e(find(QM_Rv_1e(:, 2) == GPSOS), 4);
            obs = (obs_BsRS - obs_RvRS) - (obs_BsOS - obs_RvOS);
            %% DD 계산치 생성 파트 - 기타위성 좌표 계산(기준위성 좌표는 이미 계산 완료)
            icol = PickEPH(eph, GPSOS, gs);
            STT = GetSTTbrdm(gs, eph, icol, x(1:3)); % :  OS 위성위치를 추정좌표 기준으로 변경 11/9/14
            tc = gs - STT;
            vec_GPSOS = GetSatPosNC(eph, icol, tc);
            vec_GPSOS = RotSatPos(vec_GPSOS, STT);
            %% DD 계산치 생성 파트 - 거리계산치를 각각 계산한 다음 DD 계산치 계산
            vec_BsRS = vec_GPSRS - TruePosBs_PP';  com_BsRS = norm(vec_BsRS);
            vec_RvRS = vec_GPSRS - x(1:3);    com_RvRS = norm(vec_RvRS);
            vec_BsOS = vec_GPSOS - TruePosBs_PP';  com_BsOS = norm(vec_BsOS);
            vec_RvOS = vec_GPSOS - x(1:3);    com_RvOS = norm(vec_RvOS);
            com = (com_BsRS - com_RvRS) - (com_BsOS - com_RvOS);
            y = obs*Lambda_L1 - com;
            Y(cnt2,1) = y;
            
            if Iter == 1
                %% 각 위성 az, el 저장
                [azBsRS,elBsRS] = xyz2azel(vec_BsRS, AppLatBs, AppLonBs);
                azelGPSBsRS(j,:) = [azBsRS, elBsRS];
                [azBsOS,elBsOS] = xyz2azel(vec_BsOS, AppLatBs, AppLonBs);
                azelGPSBsOS(cnt,:) = [azBsOS, elBsOS];
                %% Weighting
                SNRRS(cnt,:) = [GPS_snr_Bs_RS, GPS_snr_Rv_RS,gs,GPSRS];
                SNROS(cnt,:) = [GPS_snr_Bs_OS, GPS_snr_Rv_OS,gs,GPSOS];
                DDel(cnt,:) = [elBsRS, elBsOS,gs,GPSOS];
                W = 1;
                %             W = DDMakeW_elsnr(S1RS(cnt,:),S1OS(cnt,:),DDel(cnt,:));
                weight(cnt,:) = W;
                cnt = cnt + 1;
                el = elBsOS;
            end
            if el >= eleCut
                % H 행렬 계산 파트
%                 H(1,1) = vec_RvRS(1)/com_RvRS - vec_RvOS(1)/com_RvOS;
%                 H(1,2) = vec_RvRS(2)/com_RvRS - vec_RvOS(2)/com_RvOS;
%                 H(1,3) = vec_RvRS(3)/com_RvRS - vec_RvOS(3)/com_RvOS;
                H(cnt2,1) = vec_RvRS(1)/com_RvRS - vec_RvOS(1)/com_RvOS;
                H(cnt2,2) = vec_RvRS(2)/com_RvRS - vec_RvOS(2)/com_RvOS;
                H(cnt2,3) = vec_RvRS(3)/com_RvRS - vec_RvOS(3)/com_RvOS;
                H(cnt2,4) = 1;
                cnt2= cnt2 + 1;
%                 HTH = HTH + H'*W*H;
%                 HTy = HTy + H'*W*y;
            end
            
        end
        %% BDS DD part
        for kS = 1:length(BDS(:,1))
            if kS == BDSindexRS || BDS(kS, 3) < eleCut
                if BDS(kS, 3) < eleCut
                    NoBDSSatsUsed = NoBDSSatsUsed - 1;
                    disp([BDSRS BDS(kS, 2)])
                end
                continue
            end
            %% DD 관측치 생성 파트 --- 추후에 for 루프 밖으로 빼야 함 11/8/14
            BDSOS = BDS(kS,2);
            OtherBDSSats(kS,1) = BDSOS;
            BDS_snr_Bs_OS = QM_Bs_snr_1e(find(QM_Bs_snr_1e(:,2) == BDSOS),4);           % 베이스, BDS 기준 위성 snr 값
            BDS_snr_Rv_OS = QM_Rv_snr_1e(find(QM_Rv_snr_1e(:,2) == BDSOS),4);           % 베이스, BDS 기준 위성 snr 값
            
            S1BsOS = QM_Bs_snr_1e(find(QM_Bs_snr_1e(:,2) == BDSOS), 4);      % BASE OS SNR matrix
            S1RvOS = QM_Rv_snr_1e(find(QM_Rv_snr_1e(:,2) == BDSOS), 4);         % ROVER OS SNR matrix
            
            obs_BsRS = QM_Bs_1e(find(QM_Bs_1e(:, 2) == BDSRS), 4);
            obs_RvRS = QM_Rv_1e(find(QM_Rv_1e(:, 2) == BDSRS), 4);
            obs_BsOS = QM_Bs_1e(find(QM_Bs_1e(:, 2) == BDSOS), 4);
            obs_RvOS = QM_Rv_1e(find(QM_Rv_1e(:, 2) == BDSOS), 4);
            obs = (obs_BsRS - obs_RvRS) - (obs_BsOS - obs_RvOS);
            %% DD 계산치 생성 파트 - 기타위성 좌표 계산(기준위성 좌표는 이미 계산 완료)
            icol = PickEPH(eph, BDSOS, gs);
            STT = GetSTTbrdm(gs, eph, icol, x(1:3)); % :  OS 위성위치를 추정좌표 기준으로 변경 11/9/14
            tc = gs - STT;
            vec_BDSOS = GetSatPosNC(eph, icol, tc);
            vec_BDSOS = RotSatPos(vec_BDSOS, STT);
            %% DD 계산치 생성 파트 - 거리계산치를 각각 계산한 다음 DD 계산치 계산
            vec_BsRS = vec_BDSRS - TruePosBs_PP';  com_BsRS = norm(vec_BsRS);
            vec_RvRS = vec_BDSRS - x(1:3);    com_RvRS = norm(vec_RvRS);
            vec_BsOS = vec_BDSOS - TruePosBs_PP';  com_BsOS = norm(vec_BsOS);
            vec_RvOS = vec_BDSOS - x(1:3);    com_RvOS = norm(vec_RvOS);
            com = (com_BsRS - com_RvRS) - (com_BsOS - com_RvOS);
            y = obs -com;
            Y(cnt2,1) = y;
            if Iter == 1
                %% 각 위성 az, el 저장
                [azBsRS,elBsRS] = xyz2azel(vec_BsRS, AppLatBs, AppLonBs);
                azelGPSBsRS(j,:) = [azBsRS, elBsRS];
                [azBsOS,elBsOS] = xyz2azel(vec_BsOS, AppLatBs, AppLonBs);
                azelGPSBsOS(cnt,:) = [azBsOS, elBsOS];
                %% Weighting
                SNRRS(cnt,:) = [BDS_snr_Bs_RS, BDS_snr_Rv_RS,gs,BDSRS];
                SNROS(cnt,:) = [BDS_snr_Bs_OS, BDS_snr_Rv_OS,gs,BDSOS];
                DDel(cnt,:) = [elBsRS, elBsOS,gs,BDSOS];
                
                W =1;
                %             W = DDMakeW_elsnr(S1RS(cnt,:),S1OS(cnt,:),DDel(cnt,:));
                weight(cnt,:) = W;
                cnt = cnt + 1;
            end
            %% H 행렬 계산 파트
%             H(1,1) = vec_RvRS(1)/com_RvRS - vec_RvOS(1)/com_RvOS;
%             H(1,2) = vec_RvRS(2)/com_RvRS - vec_RvOS(2)/com_RvOS;
%             H(1,3) = vec_RvRS(3)/com_RvRS - vec_RvOS(3)/com_RvOS;
            H(cnt2,1) = vec_RvRS(1)/com_RvRS - vec_RvOS(1)/com_RvOS;
            H(cnt2,2) = vec_RvRS(2)/com_RvRS - vec_RvOS(2)/com_RvOS;
            H(cnt2,3) = vec_RvRS(3)/com_RvRS - vec_RvOS(3)/com_RvOS;
            H(cnt2,4) = 1;
            cnt2= cnt2 + 1;
%             HTH = HTH + H'*W*H;
%             HTy = HTy + H'*W*y;
        end
        
%         OTHER{j,1} = OtherSats;
%         OTHER{j,2} = indxUsedSat;
            HTH = H'*W*H;
            HTy = H'*W*Y;

        xhat = inv(HTH) * HTy;
        x = x + xhat;
        XHAT(j,1:2) =[gs, norm(xhat)];
        if norm(xhat(1:3)) < EpsStop;
            nEst = nEst + 1;
            estm(nEst,1) =gs;
            estm(nEst,2:4) =x(1:3);
            estm(nEst,5) = length(GPS(:,1));        % GPS 관측 위성 수
            estm(nEst,6) = length(BDS(:,1));        % BDS 관측 위성 수
            estm(nEst,7) = NoGPSSatsUsed;           % : snr issue 가 없을때
            estm(nEst,8) = NoBDSSatsUsed;           % : snr issue 가 없을때
            clear H;
            clear HTH;
            clear HTy;
            clear Y;
            
            
            break;
        end
        
    end
   
    %% rover's longi, lati
    rover_gd(nEst,1) = gs;
    rover_gd(nEst,2:4) = xyz2gd(estm(nEst,2:4)); % rover의 xyz값 gd 로 변환
    AppLat = rover_gd(nEst,2); AppLon = rover_gd(nEst,3);
    
end

%% 측위오차 분석 & 그래프 작성
% estm = estm(1:nEst, :);

[DDdXYZ, DDdNEV, DDdis, DDrms, DDstd] = PostErrorsDDmv2(estm, Base, Truedis,8, 12, visiSat);    % 임의의 장소에서 이동 측위시

%% 각 위성별 SNR Plot
% DDPlotQM(QM11, QM22, 141, 'Base', 'Rover')








