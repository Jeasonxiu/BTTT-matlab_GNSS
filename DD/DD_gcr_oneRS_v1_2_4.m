%% 코드의사거리 이중차분 알고리즘
% 07/01/2016 : Joonseong
close all; clear all;
% % load('DD_gc_test.mat');
% 
% %% 항법 시스템 선택 (1=GPS, 2=BDS, 3=GLO, 4=GPS/BDS, 5=GPS/GLO, 6=GLO/BDS, 7=GPS/GLO/BDS)
SYS = 7;
% %% Static or Kinematic 선택 (1=kinematic, 0=static)
% Dynamic = 0;
% if Dynamic == 0
%     Truedis = 9.9209;           % 대성 디폴리스 A,B
% elseif Dynamic == 1
%     Base_vrs = load('DD1Bs_adm.txt');
%     Rover_vrs = load('DD1Rv_adm.txt');
%     %     Base_vrs = load('DD2Bs_adm.txt');
%     %     Rover_vrs = load('DD2Rv_adm.txt');
%     %     Base_vrs = load('PTCO4_joon_170216_adm.txt');
%     %     Rover_vrs = load('PTCO4_hyunu_170216_adm.txt');
% end
% 
% %% 계산 날짜 설정
% DOY = 025; YY  = 17;        % 대성 AB
% % DOY = 047; YY  = 17;
% % DOY = 081; YY  = 17;        % 오송
% [gw, gd] = ydoy2gwgd(YY, DOY); %: GPS WEEK 결정
% 
% %% QM 파일 핸들링
% QMfileBs = 'QM170125_A';        % 대성 A
% QMfileRv = 'QM170125_B';        % 대성 B
% % QMfileBs = 'QMfile_Bs_8';
% % QMfileRv = 'QMfile_Rv_8';
% % QMfileBs = 'ublox_joon';
% % QMfileRv = 'ublox_hyunu';
% % QMfileBs = 'QM170322_Bs_1';        % 오송 bs 1
% % QMfileRv = 'QM170322_Rv_1';        % 오송 rv 1
% % QMfileBs = 'QM170322_Bs_2';        % 오송 bs 2
% % QMfileRv = 'QM170322_Rv_2';        % 오송 rv 2
% 
% %% 좌표 참값
% TruePos = [-3041235.57800000,4053941.67700000,3859881.01300000];        % 대성 A
% % TruePos = [-3041241.74100000,4053944.14300000,3859873.64000000];        % 대성 B
% % TruePos = [-3027386.463213997 4071581.638074351 3852036.292033684]; % 현우집앞
% % TruePos = [-3108706.97103816,4078522.84147406,3779757.23816543];    % 오송 1
% % TruePos = [-3108697.15866998,4078501.37490046,3779789.12574991];    % 오송 2
% %% Bs와 EPH를 불러올 경우
% Bs = [];
% % load('DD_gc_Bs.mat');         % 대성 AB
% % load('DD_gc_Bs_joon_m.mat');
% % load('DD_170322_gc_Bs_1.mat');
% % load('DD_170322_gc_Bs_2.mat');
% 
% %% 불변 변수 설정 : 빛의 속도, 관측치
% CCC = 299792458.;   % CCC = Speed of Light [m/s]
% Freq_L1 = 1575.42e6;        % L1 주파수
% Lambda_L1 = CCC/Freq_L1;    % L1 파장
% ObsType = 120;      % 사용할 관측지 설정 - 120 : C/A = C1
% ObsType2 = 141;      % 사용할 관측치 설정 - 141: snr = S1
% %% 임계고도각 설정
% eleCut = 15;
% %% 윤초 핸들링
% LeapSec = 18;
% LeapSecBDS = 14;
% %% True Distance
% % Truedis = 1.41;         % 신항대로 실험
% 
% %% 항법메시지 호출
% navfile = strcat('brdm', num2str(DOY,'%03d'), '0.', num2str(YY,'%02d'), 'p');   %: Navigation RINEX file
% %
% TauC = ReadTauC2(navfile);
% 
% %% 사용할 관측치 선정
% g_ObsType = 120; % gps C1
% g_ObsType_snr = 141;
% g_ObsType_L1 = 111;    % L1
% c_ObsType = 220; % bds C1
% c_ObsType_snr = 241;
% c_ObsType_L1 = 211;    % L1
% r_ObsType = 320; % bds C1
% r_ObsType_snr = 341;
% r_ObsType_L1 = 311;    % L1
% 
% % load('modifiedQM.mat');
% %% QM 파일 읽어들여서 행렬로 저장하고, 사용할 관측치 추출
% %% 베이스 QM
% [arrQM_Bs, FinalPRNs_Bs, FinalTTs_Bs] = ReadQM(QMfileBs);
% % arrQM_Bs = arrQM_Bs(find(arrQM_Bs(:,3) < 300),:);
% arrQM_Bs(:,1) = round(arrQM_Bs(:,1));
% % QM_Bs = SelectQM_gc(arrQM_Bs, g_ObsType, c_ObsType);
% % QM_Bs_snr = SelectQM_gc(arrQM_Bs, g_ObsType_snr, c_ObsType_snr);
% % QM_Bs_L1 = SelectQM_gc(arrQM_Bs, g_ObsType_L1, c_ObsType_L1);
% QM_Bs = SelectQM_gcr(arrQM_Bs, g_ObsType, c_ObsType, r_ObsType);
% QM_Bs_snr = SelectQM_gcr(arrQM_Bs, g_ObsType_snr, c_ObsType_snr, r_ObsType_snr);
% QM_Bs_L1 = SelectQM_gcr(arrQM_Bs, g_ObsType_L1, c_ObsType_L1, r_ObsType_L1);
% cp_prn_Bs = unique(QM_Bs(:,2));
% 
% %% 로버 QM
% [arrQM_Rv, FinalPRNs_Rv, FinalTTs_Rv] = ReadQM(QMfileRv);
% % arrQM_Rv = arrQM_Rv(find(arrQM_Rv(:,3) < 300),:);
% arrQM_Rv(:,1) = round(arrQM_Rv(:,1));
% % QM_Rv = SelectQM_gc(arrQM_Rv, g_ObsType, c_ObsType);
% % QM_Rv_snr = SelectQM_gc(arrQM_Rv, g_ObsType_snr, c_ObsType_snr);
% % QM_Rv_L1 = SelectQM_gc(arrQM_Rv, g_ObsType_L1, c_ObsType_L1);
% QM_Rv = SelectQM_gcr(arrQM_Rv, g_ObsType, c_ObsType, r_ObsType);
% QM_Rv_snr = SelectQM_gcr(arrQM_Rv, g_ObsType_snr, c_ObsType_snr, r_ObsType_snr);
% QM_Rv_L1 = SelectQM_gcr(arrQM_Rv, g_ObsType_L1, c_ObsType_L1, r_ObsType_L1);
% cp_prn_Rv = unique(QM_Rv(:,2));
% 
% %% 항법메시지를 읽어들여서 행렬로 저장하고, Klobuchar 모델 추출
% gps_nav = strcat(navfile(1:3),'c',navfile(5:8),'.',navfile(10:11),'n');
% fid = fopen(gps_nav,'r');
% if fid == -1
%     al = zeros(4,1); be = zeros(4,1);
% else
%     [al, be] = GetALBE(gps_nav);
% end
% 
% %% 라이넥스 파일에서 대략적인 관측소 좌표를 뽑아냄
% AppPos = TruePos;
% 
% %% 라이넥스 파일에서 Base Station의 좌표를 계산함
% if isempty(Bs)
%     [eph, trashPRN, trashT]=ReadEPH_all(navfile);
%     if SYS == 1
%         disp('GPS only')
%         Bs = PP_g(QMfileBs,eph,TruePos,DOY,YY);              % with Correction
%     elseif SYS == 2
%         disp('BDS only')
%         Bs = PP_c(QMfileBs,eph,TruePos,DOY,YY);              % with Correction
%     elseif SYS == 3
%         disp('GLO only')
%         Bs = PP_r(QMfileBs,eph,TruePos,DOY,YY);              % with Correction
%     elseif SYS == 4
%         disp('GPS BDS')
%         Bs = PP_gc(QMfileBs,eph,TruePos,DOY,YY);              % with Correction
%         %         Bs = PP_gc_v1(QMfileBs,eph,TruePos,DOY,YY,0);              % with Correction
%     elseif SYS == 5
%         disp('GPS GLO')
%         Bs = PP_gr(QMfileBs,eph,TruePos,DOY,YY);              % with Correction
%     elseif SYS == 6
%         disp('BDS GLO')
%         Bs = PP_rc(QMfileBs,eph,TruePos,DOY,YY);              % with Correction
%     elseif SYS == 7
%         disp('GPS BDS GLO')
%         Bs = PP_grc(QMfileBs,eph,TruePos,DOY,YY);              % with Correction
%     end
% end
% ephGLO = eph;
% ephGLO=ephGLO(ephGLO(:,1)<400&ephGLO(:,1)>300,:);ephGLO(:,1)=ephGLO(:,1)-300;
% Bs(:,1) = round(Bs(:,1));
% %% 두 QM 파일에서 공통시각(epoch) 추출
% FinalTTs = intersect(Bs(:, 1), QM_Rv(:, 1));
% 
% 
% %% 추정에 필요한 초기치 설정
% MaxIter = 10;
% EpsStop = 1e-4;
% x = AppPos';
% if SYS == 4 | SYS == 5 | SYS == 6
%     x = [AppPos' ; 0];
% elseif SYS ==7
%     x = [AppPos' ; 0; 0];
% end
% 
% deltat = 5;
% 
% %% 추정과정 시작
% NoEpochs = length(FinalTTs);
% % estm = zeros(NoEpochs, 4);  % : c1(gs), c2/3/4(X/Y/Z of Rv), c5(#Sats_obs), c6(#Sats_used)
% SatPosArr_before=0;icolArr_before=zeros(24,2);
% nEst = 0;
switch SYS
    case 1
        load('DD_PP_g_170125_wc.mat')
    case 2
        load('DD_PP_c_170125_wc.mat')
    case 3
        load('DD_PP_r_170125_wc.mat')
    case 4
        load('DD_PP_gc_170125_wc.mat')
    case 5
        load('DD_PP_gr_170125_wc.mat')
    case 6
        load('DD_PP_rc_170125_wc.mat')
    case 7
        load('DD_PP_grc_170125_wc.mat')
        MaxIter = 4;
        eleCut = 25;
end
cnt = 1;

for j = 1:NoEpochs
    % for j = 100:500
    j
    gs = FinalTTs(j);
    %% 해당 시각 Bs의 위치 값 찾기(PP-LS)
    Base(j,:) = Bs(find(Bs(:,1) == gs),:);
    TruePosBs_PP = Base(j,2:4);                                             % 베이스 좌표를 기준좌표로 설정
    %         TruePosBs_PP = TruePos;                                             % 고정 좌표를 기준좌표로 설정
    %% base의 gd
    base_gd(j,1) = gs;
    base_gd(j,2:4) = xyz2gd(TruePosBs_PP(1:3));                             % Base의 xyz값 gd 로 변환
    rover_gd(j,1:4) = [gs, xyz2gd(x(1:3)')];
    LatBs = base_gd(j,2); LonBs = base_gd(j,3);
    LatRv = rover_gd(j,2); LonRv = rover_gd(j,3);
    %% 해당 시각 gs의 베이스 관측치 선별
    indexQM1 = find(QM_Bs(:,1) == gs);
    QM_Bs_1e = QM_Bs(indexQM1,:);                                           % Base Pseudo-Range(gs)
    indexQM1 = find(QM_Bs_snr(:,1) == gs);
    QM_Bs_snr_1e = QM_Bs_snr(indexQM1,:);                                   % Base SNR(gs)
    indexQM1 = find(QM_Bs_L1(:,1) == gs);
    QM_Bs_L1_1e = QM_Bs_L1(indexQM1,:);                                   % Base L1(gs)
    
    %% 해당 시각 gs의 로버 관측치 선별
    indexQM2 = find(QM_Rv(:,1) == gs);
    QM_Rv_1e = QM_Rv(indexQM2,:);                                           % Rover Pseudo-Range(gs)
    indexQM2 = find(QM_Rv_snr(:,1) == gs);
    QM_Rv_snr_1e = QM_Rv_snr(indexQM2,:);                                   % Rover SNR(gs)
    indexQM2 = find(QM_Rv_L1(:,1) == gs);
    QM_Rv_L1_1e = QM_Rv_L1(indexQM2,:);                                   % Base L1(gs)
    %% 선택된 system 별, 한 epoch의 QM 행렬 선택
    switch SYS
        case 1
            QM_Bs_1e = QM_Bs_1e(find(QM_Bs_1e(:,2) < 200),:);
            QM_Rv_1e = QM_Rv_1e(find(QM_Rv_1e(:,2) < 200),:);
            QM_Bs_snr_1e = QM_Bs_snr_1e(find(QM_Bs_snr_1e(:,2) < 200),:);
            QM_Rv_snr_1e = QM_Rv_snr_1e(find(QM_Rv_snr_1e(:,2) < 200),:);
            QM_Bs_L1_1e = QM_Bs_L1_1e(find(QM_Bs_L1_1e(:,2) < 200),:);
            QM_Rv_L1_1e = QM_Rv_L1_1e(find(QM_Rv_L1_1e(:,2) < 200),:);
            Sats = intersect(QM_Bs_1e(:, 2), QM_Rv_1e(:, 2));
            Sats = intersect(Sats, unique(eph(find(eph(:,22) == 0),18)));
            NoSats = length(Sats);
            if NoSats > 5
                %% 기준위성 RS와 다른위성 OS 선별/ SatsEl - c1(gs), c2(prn), c3(el)
                [SatsEl, indxRS] = PickRSel_grc(gs, Sats, eph, TruePosBs_PP);  % : RS 기준위성
            end
        case 2
            QM_Bs_1e = QM_Bs_1e(find(QM_Bs_1e(:,2) > 200 & QM_Bs_1e(:,2) < 300),:);
            QM_Rv_1e = QM_Rv_1e(find(QM_Rv_1e(:,2) > 200 & QM_Rv_1e(:,2) < 300),:);
            QM_Bs_snr_1e = QM_Bs_snr_1e(find(QM_Bs_snr_1e(:,2) > 200 & QM_Bs_snr_1e(:,2) < 300),:);
            QM_Rv_snr_1e = QM_Rv_snr_1e(find(QM_Rv_snr_1e(:,2) > 200 & QM_Rv_snr_1e(:,2) < 300),:);
            QM_Bs_L1_1e = QM_Bs_L1_1e(find(QM_Bs_L1_1e(:,2) > 200 & QM_Bs_L1_1e(:,2) < 300),:);
            QM_Rv_L1_1e = QM_Rv_L1_1e(find(QM_Rv_L1_1e(:,2) > 200 & QM_Rv_L1_1e(:,2) < 300),:);
            Sats = intersect(QM_Bs_1e(:, 2), QM_Rv_1e(:, 2));
            Sats = intersect(Sats, unique(eph(find(eph(:,22) == 0),18)));
            NoSats = length(Sats);
            if NoSats > 5
                %% 기준위성 RS와 다른위성 OS 선별/ SatsEl - c1(gs), c2(prn), c3(el)
                [SatsEl, indxRS] = PickRSel_grc(gs, Sats, eph, TruePosBs_PP);  % : RS 기준위성
            end
        case 3
            QM_Bs_1e = QM_Bs_1e(find(QM_Bs_1e(:,2) > 300),:);
            QM_Rv_1e = QM_Rv_1e(find(QM_Rv_1e(:,2) > 300),:);
            QM_Bs_snr_1e = QM_Bs_snr_1e(find(QM_Bs_snr_1e(:,2) > 300),:);
            QM_Rv_snr_1e = QM_Rv_snr_1e(find(QM_Rv_snr_1e(:,2) > 300),:);
            QM_Bs_L1_1e = QM_Bs_L1_1e(find(QM_Bs_L1_1e(:,2) > 300),:);
            QM_Rv_L1_1e = QM_Rv_L1_1e(find(QM_Rv_L1_1e(:,2) > 300),:);
            Sats = intersect(QM_Bs_1e(:, 2), QM_Rv_1e(:, 2));
            Sats = intersect(Sats, unique(eph(find(eph(:,22) == 0),18)));
            NoSats = length(Sats);
            if NoSats > 5
                %% 기준위성 RS와 다른위성 OS 선별/ SatsEl - c1(gs), c2(prn), c3(el)
                [SatsEl, indxRS] = PickRSel_grc(gs, Sats, eph, TruePosBs_PP);  % : RS 기준위성
            end
        case 4
            QM_Bs_1e = QM_Bs_1e(find(QM_Bs_1e(:,2) < 300),:);
            QM_Rv_1e = QM_Rv_1e(find(QM_Rv_1e(:,2) < 300),:);
            QM_Bs_snr_1e = QM_Bs_snr_1e(find(QM_Bs_snr_1e(:,2) < 300),:);
            QM_Rv_snr_1e = QM_Rv_snr_1e(find(QM_Rv_snr_1e(:,2) < 300),:);
            QM_Bs_L1_1e = QM_Bs_L1_1e(find(QM_Bs_L1_1e(:,2) < 300),:);
            QM_Rv_L1_1e = QM_Rv_L1_1e(find(QM_Rv_L1_1e(:,2) < 300),:);
            Sats = intersect(QM_Bs_1e(:, 2), QM_Rv_1e(:, 2));
            Sats = intersect(Sats, unique(eph(find(eph(:,22) == 0),18)));
            NoSats = length(Sats);
            if NoSats >= 6
                %% 기준위성 RS와 다른위성 OS 선별/ SatsEl - c1(gs), c2(prn), c3(el)
                [SatsEl, indxRS] = PickRSel_grc(gs, Sats, eph, TruePosBs_PP);  % : RS 기준위성
            end
        case 5
            QM_Bs_1e = QM_Bs_1e(find(QM_Bs_1e(:,2) < 200 & QM_Bs_1e(:,2) > 300),:);
            QM_Rv_1e = QM_Rv_1e(find(QM_Rv_1e(:,2) < 200 & QM_Rv_1e(:,2) > 300),:);
            QM_Bs_snr_1e = QM_Bs_snr_1e(find(QM_Bs_snr_1e(:,2) < 200 & QM_Bs_snr_1e(:,2) > 300),:);
            QM_Rv_snr_1e = QM_Rv_snr_1e(find(QM_Rv_snr_1e(:,2) < 200 & QM_Rv_snr_1e(:,2) > 300),:);
            QM_Bs_L1_1e = QM_Bs_L1_1e(find(QM_Bs_L1_1e(:,2) < 200 & QM_Bs_L1_1e(:,2) > 300),:);
            QM_Rv_L1_1e = QM_Rv_L1_1e(find(QM_Rv_L1_1e(:,2) < 200 & QM_Rv_L1_1e(:,2) > 300),:);
            Sats = intersect(QM_Bs_1e(:, 2), QM_Rv_1e(:, 2));
            Sats = intersect(Sats, unique(eph(find(eph(:,22) == 0),18)));
            NoSats = length(Sats);
            if NoSats > 5
                %% 기준위성 RS와 다른위성 OS 선별/ SatsEl - c1(gs), c2(prn), c3(el)
                [SatsEl, indxRS] = PickRSel_grc(gs, Sats, eph, TruePosBs_PP);  % : RS 기준위성
            end
        case 6
            QM_Bs_1e = QM_Bs_1e(find(QM_Bs_1e(:,2) > 200),:);
            QM_Rv_1e = QM_Rv_1e(find(QM_Rv_1e(:,2) > 200),:);
            QM_Bs_snr_1e = QM_Bs_snr_1e(find(QM_Bs_snr_1e(:,2) > 200),:);
            QM_Rv_snr_1e = QM_Rv_snr_1e(find(QM_Rv_snr_1e(:,2) > 200),:);
            QM_Bs_L1_1e = QM_Bs_L1_1e(find(QM_Bs_L1_1e(:,2) > 200),:);
            QM_Rv_L1_1e = QM_Rv_L1_1e(find(QM_Rv_L1_1e(:,2) > 200),:);
            Sats = intersect(QM_Bs_1e(:, 2), QM_Rv_1e(:, 2));
            Sats = intersect(Sats, unique(eph(find(eph(:,22) == 0),18)));
            NoSats = length(Sats);
            if NoSats >= 6
                %% 기준위성 RS와 다른위성 OS 선별/ SatsEl - c1(gs), c2(prn), c3(el)
                [SatsEl, indxRS] = PickRSel_grc(gs, Sats, eph, TruePosBs_PP);  % : RS 기준위성
            end
        case 7
            Sats = intersect(QM_Bs_1e(:, 2), QM_Rv_1e(:, 2));
            Sats = intersect(Sats, unique(eph(find(eph(:,22) == 0),18)));
            NoSats = length(Sats);
            if NoSats >= 6
                %% 기준위성 RS와 다른위성 OS 선별/ SatsEl - c1(gs), c2(prn), c3(el)
                [SatsEl, indxRS] = PickRSel_grc(gs, Sats, eph, TruePosBs_PP);  % : RS 기준위성
            end
    end
    if NoSats >= 6
        
        GPS = SatsEl(find(SatsEl(:,2) < 200),:);                                % SatsEl 중 GPS 위성군 선택
        BDS = SatsEl(find(SatsEl(:,2) > 200 & SatsEl(:,2) < 300),:);            % SatsEl 중 BDS 위성군 선택
        GLO = SatsEl(find(SatsEl(:,2) > 300),:);                                % SatsEl 중 BDS 위성군 선택
        RS = Sats(indxRS); RefSV(j,1) = RS;
        RefSV(j,1:2) = [gs, RS];                               % 기준위성 PRN 저장
        
        %% gs, 공통위성 수, Base GPS, Rover GPS, Base GLO, Rover GLO...
        %% [gs, 전체 관측 위성수, 베이스 GPS 관측 수, 로버 GPS 관측수,...
        %%  베이스 BDS 관측 수, 로버 BDS 관측 수, 베이스 전체 관측수, 로버 전체 관측수]
        visiSat(j,1) = gs; visiSat(j,2) = length(GPS(:,1)) + length(BDS(:,1)) + length(GLO(:,1));
        visiSat(j,3) = length(QM_Bs_1e(find(QM_Bs_1e(:,2) < 200), 2));
        visiSat(j,4) = length(QM_Rv_1e(find(QM_Rv_1e(:,2) < 200), 2));
        visiSat(j,5) = length(QM_Bs_1e(find(QM_Bs_1e(:,2) > 200 & QM_Bs_1e(:,2) < 300), 2));
        visiSat(j,6) = length(QM_Rv_1e(find(QM_Rv_1e(:,2) > 200 & QM_Rv_1e(:,2) < 300), 2));
        visiSat(j,7) = length(QM_Bs_1e(find(QM_Bs_1e(:,2) > 300), 2));
        visiSat(j,8) = length(QM_Rv_1e(find(QM_Rv_1e(:,2) > 300), 2));
        visiSat(j,9) = visiSat(j,3)+visiSat(j,5)+visiSat(j,7);
        visiSat(j,10) = visiSat(j,4)+visiSat(j,6)+visiSat(j,8);
        
        %% 기준 위성 snr 값 추출
        snr_Bs_RS = QM_Bs_snr_1e(find(QM_Bs_snr_1e(:,2) == RS),4);              % 베이스, 기준 위성 snr 값
        snr_Rv_RS = QM_Rv_snr_1e(find(QM_Rv_snr_1e(:,2) == RS),4);              % 로버, 기준 위성 snr 값
        L1_Bs_RS = QM_Bs_L1_1e(find(QM_Bs_L1_1e(:,2) == RS),4);                 % 베이스, GPS 기준 위성 snr 값
        L1_Rv_RS = QM_Rv_L1_1e(find(QM_Rv_L1_1e(:,2) == RS),4);                 % 로버, GPS 기준 위성 snr 값
        %% Iteration 시작
        for Iter = 1:MaxIter
            cnt2=1;
            NoGPSSatsUsed = length(GPS(:,1));                                   % GPS 사용 위성수 초기값
            NoBDSSatsUsed = length(BDS(:,1));                                   % BDS 사용 위성수 초기값
            NoGLOSatsUsed = length(GLO(:,1));                                   % GLO 사용 위성수 초기값
            NoSatsUsed = length(Sats(:,1));                                     % 전체 사용 위성수 초기값
            
            %% 선택 시스템 추정
            if SYS == 1 | SYS == 2
                %% 기준위성 좌표 먼저 계산 - Bs 기준
                icol = PickEPH(eph, RS ,gs);
                t_delay = QM_Bs_1e(find(QM_Bs_1e(:,2) == RS),4)/CCC;
                t_delay = (QM_Rv_1e(find(QM_Rv_1e(:,2) == RS),4)/(CCC+t_delay))/2;
                tc = gs - t_delay;
                vec_RS = GetSatPosNC_GC(eph, icol, tc);
                vec_RS = RotSatPos(vec_RS, t_delay);                                    % 지구 자전효과 고려
                %% GPS, BDS or GPS/BDS DD part
                for kS = 1:NoSats
                    if kS == indxRS || SatsEl(kS, 3) < eleCut
                        if SatsEl(kS, 3) < eleCut
                            NoSatsUsed = NoSatsUsed - 1;
                            if SatsEl(kS, 3) < 200
                                NoGPSSatsUsed = NoGPSSatsUsed - 1;
                            elseif SatsEl(kS, 3) > 200 & SatsEl(kS, 3) < 300
                                NoBDSSatsUsed = NoBDSSatsUsed - 1;
                            elseif SatsEl(kS, 3) > 300
                                NoGLOSatsUsed = NoGLOSatsUsed - 1;
                            end
                            disp([RS SatsEl(kS, 2)])
                        end
                        continue
                    end
                    %% DD 관측치 생성 파트
                    OS = Sats(kS,1);
                    OtherSats(kS,1) = OS;
                    snr_Bs_OS = QM_Bs_snr_1e(find(QM_Bs_snr_1e(:,2) == OS),4);           % 베이스, GPS 기준 위성 snr 값
                    snr_Rv_OS = QM_Rv_snr_1e(find(QM_Rv_snr_1e(:,2) == OS),4);           % 베이스, GPS 기준 위성 snr 값
                    
                    L1_Bs_OS = QM_Bs_L1_1e(find(QM_Bs_L1_1e(:,2) == OS),4);           % 베이스, GPS 기준 위성 snr 값
                    L1_Rv_OS = QM_Rv_L1_1e(find(QM_Rv_L1_1e(:,2) == OS),4);           % 베이스, GPS 기준 위성 snr 값
                    
                    obs_BsRS = QM_Bs_1e(find(QM_Bs_1e(:, 2) == RS), 4);
                    obs_RvRS = QM_Rv_1e(find(QM_Rv_1e(:, 2) == RS), 4);
                    obs_BsOS = QM_Bs_1e(find(QM_Bs_1e(:, 2) == OS), 4);
                    obs_RvOS = QM_Rv_1e(find(QM_Rv_1e(:, 2) == OS), 4);
                    
                    obs = (obs_BsRS - obs_RvRS) - (obs_BsOS - obs_RvOS);
                    obs2 = (L1_Bs_RS - L1_Rv_RS) - (L1_Bs_OS - L1_Rv_OS);
                    
                    OBS(cnt,1:2) = [(obs_BsRS - obs_RvRS) (obs_BsOS - obs_RvOS)];
                    %                         OBS2(cnt,1:2) = [(L1_Bs_RS - L1_Rv_RS) (L1_Bs_OS - L1_Rv_OS)];
                    %% DD 계산치 생성 파트 - 기타위성 좌표 계산(기준위성 좌표는 이미 계산 완료)
                    icol = PickEPH(eph, OS, gs);
                    t_delay = QM_Bs_1e(find(QM_Bs_1e(:,2) == OS),4)/CCC;
                    t_delay = (QM_Rv_1e(find(QM_Rv_1e(:,2) == OS),4)/(CCC+t_delay))/2;
                    tc = gs - t_delay;
                    vec_OS = GetSatPosNC_GC(eph, icol, tc);
                    vec_OS = RotSatPos(vec_OS, t_delay);
                    %% DD 계산치 생성 파트 - 거리계산치를 각각 계산한 다음 DD 계산치 계산
                    vec_BsRS = vec_RS - TruePosBs_PP';  com_BsRS = norm(vec_BsRS);
                    vec_RvRS = vec_RS - x(1:3);    com_RvRS = norm(vec_RvRS);
                    vec_BsOS = vec_OS - TruePosBs_PP';  com_BsOS = norm(vec_BsOS);
                    vec_RvOS = vec_OS - x(1:3);    com_RvOS = norm(vec_RvOS);
                    %                     vec_BsRS = TruePosBs_PP' - vec_RS;  com_BsRS = norm(vec_BsRS);
                    %                     vec_RvRS = x(1:3) - vec_RS;    com_RvRS = norm(vec_RvRS);
                    %                     vec_BsOS = TruePosBs_PP' - vec_OS;  com_BsOS = norm(vec_BsOS);
                    %                     vec_RvOS = x(1:3) - vec_OS;    com_RvOS = norm(vec_RvOS);
                    com = (com_BsRS - com_RvRS) - (com_BsOS - com_RvOS);
                    y = obs - com;
                    Y(cnt2,1) = y;
                    
                    if Iter == 1
                        %% 각 위성 az, el 저장
                        [azBsRS,elBsRS] = xyz2azel(vec_BsRS, LatBs, LonBs);
                        azelBsRS(j,:) = [azBsRS, elBsRS, RS];
                        [azBsOS,elBsOS] = xyz2azel(vec_BsOS, LatBs, LonBs);
                        [azRvOS,elRvOS] = xyz2azel(vec_RvOS, LatRv, LonRv);
                        azelBsOS(cnt,:) = [azBsOS, elBsOS, OS];
                        %% Weighting
                        SNRRS(j,:) = [snr_Bs_RS, snr_Rv_RS,gs,RS];
                        SNROS(cnt,:) = [snr_Bs_OS, snr_Rv_OS,gs,OS];
                        DDel(cnt,:) = [elBsRS, elBsOS,gs,OS];
                        cnt = cnt + 1;
                        el = elBsOS;
                    end
                    if el >= eleCut
                        % 가중치 계산 파트
                        RS_snr = [snr_Bs_RS, snr_Rv_RS,gs,RS];
                        OS_snr = [snr_Bs_OS, snr_Rv_OS,gs,OS];
                        if cnt2 == 1
                            elsnr_Bs(cnt2,1:2) = [elBsRS, snr_Bs_RS];
                            elsnr_Rv(cnt2,1:2) = [elBsRS, snr_Rv_OS];
                        end
                        elsnr_Bs(cnt2+1,1:2) = [elBsOS, snr_Bs_OS];
                        elsnr_Rv(cnt2+1,1:2) = [elRvOS, snr_Rv_OS];
                        H(cnt2,1) = vec_RvRS(1)/com_RvRS - vec_RvOS(1)/com_RvOS;
                        H(cnt2,2) = vec_RvRS(2)/com_RvRS - vec_RvOS(2)/com_RvOS;
                        H(cnt2,3) = vec_RvRS(3)/com_RvRS - vec_RvOS(3)/com_RvOS;
                        
                        cnt2= cnt2 + 1;
                    end
                    
                end
            elseif SYS == 3
                icol = PickEPH_GLO2(ephGLO, RS-300, gs);
                t_delay = QM_Bs_1e(find(QM_Bs_1e(:,2) == RS),4)/CCC;
                t_delay = (QM_Rv_1e(find(QM_Rv_1e(:,2) == RS),4)/(CCC+t_delay))/2;
                tc = gs - t_delay;
                tc = tc - LeapSec - TauC;
                [vec_RS,SatVel,SatLS] = GetSatPosGLO(ephGLO,icol,tc,deltat); % 방송궤도력으로 위성위치 계산
                vec_RS = RotSatPos(vec_RS, t_delay); % 지구자전효과 고려
                
                %% GLO DD part
                for kS = 1:NoSats
                    if kS == indxRS || SatsEl(kS, 3) < eleCut
                        if SatsEl(kS, 3) < eleCut
                            NoSatsUsed = NoSatsUsed - 1;
                            if SatsEl(kS, 3) < 200
                                NoGPSSatsUsed = NoGPSSatsUsed - 1;
                            elseif SatsEl(kS, 3) > 200 & SatsEl(kS, 3) < 300
                                NoBDSSatsUsed = NoBDSSatsUsed - 1;
                            elseif SatsEl(kS, 3) > 300
                                NoGLOSatsUsed = NoGLOSatsUsed - 1;
                            end
                            disp([RS SatsEl(kS, 2)])
                        end
                        continue
                    end
                    %% DD 관측치 생성 파트
                    OS = Sats(kS,1);
                    OtherSats(kS,1) = OS;
                    snr_Bs_OS = QM_Bs_snr_1e(find(QM_Bs_snr_1e(:,2) == OS),4);           % 베이스, GPS 기준 위성 snr 값
                    snr_Rv_OS = QM_Rv_snr_1e(find(QM_Rv_snr_1e(:,2) == OS),4);           % 베이스, GPS 기준 위성 snr 값
                    
                    L1_Bs_OS = QM_Bs_L1_1e(find(QM_Bs_L1_1e(:,2) == OS),4);           % 베이스, GPS 기준 위성 snr 값
                    L1_Rv_OS = QM_Rv_L1_1e(find(QM_Rv_L1_1e(:,2) == OS),4);           % 베이스, GPS 기준 위성 snr 값
                    
                    obs_BsRS = QM_Bs_1e(find(QM_Bs_1e(:, 2) == RS), 4);
                    obs_RvRS = QM_Rv_1e(find(QM_Rv_1e(:, 2) == RS), 4);
                    obs_BsOS = QM_Bs_1e(find(QM_Bs_1e(:, 2) == OS), 4);
                    obs_RvOS = QM_Rv_1e(find(QM_Rv_1e(:, 2) == OS), 4);
                    
                    obs = (obs_BsRS - obs_RvRS) - (obs_BsOS - obs_RvOS);
                    obs2 = (L1_Bs_RS - L1_Rv_RS) - (L1_Bs_OS - L1_Rv_OS);
                    
                    OBS(cnt,1:2) = [(obs_BsRS - obs_RvRS) (obs_BsOS - obs_RvOS)];
                    %                         OBS2(cnt,1:2) = [(L1_Bs_RS - L1_Rv_RS) (L1_Bs_OS - L1_Rv_OS)];
                    %% DD 계산치 생성 파트 - 기타위성 좌표 계산(기준위성 좌표는 이미 계산 완료)
                    icol = PickEPH_GLO2(ephGLO, OS-300, gs);
                    icolArr(kS,1)=OS; icolArr(kS,2)=icol;
                    TauN = ephGLO(icol,12); % : tau & gamma - 시계오차 보정에 필요
                    GammaN = ephGLO(icol,13);
                    ch_num = ephGLO(icol,16); % : channel number - 전리층 보정에 필요
                    jcol = find(SatPosArr_before(:,1)==OS);
                    if isempty(jcol); jcol=0;end
                    jcol=jcol(1);
                    %% 위성위치 계산 파트
                    if icol == icolArr_before(icolArr_before(:,1)==OS,2) % icol 변화 없으면 이전 epoch으로 계산
                        %                 disp(1)
                        t_delay = QM_Bs_1e(find(QM_Bs_1e(:,2) == OS),4)/CCC;
                        t_delay = (QM_Rv_1e(find(QM_Rv_1e(:,2) == OS),4)/(CCC+t_delay))/2;
                        tc = gs - t_delay;
                        tc = tc - LeapSec - TauC;
                        [SatPos,SatVel,SatLS] = GetSatPosGLO_ver3(SatPosArr_before,jcol,tc,deltat); % 이전 epoch값으로 위성위치 계산
                        vec_OS = RotSatPos(SatPos, t_delay); % 지구자전효과 고려
                    else % icol 변화시(새위성 또는 새로운 방송궤도력) 방송궤도력으로 계산
                        %                 disp(2)
                        t_delay = QM_Bs_1e(find(QM_Bs_1e(:,2) == OS),4)/CCC;
                        t_delay = (QM_Rv_1e(find(QM_Rv_1e(:,2) == OS),4)/(CCC+t_delay))/2;
                        tc = gs - t_delay;
                        tc = tc - LeapSec - TauC;
                        [SatPos,SatVel,SatLS] = GetSatPosGLO(ephGLO,icol,tc,deltat); % 방송궤도력으로 위성위치 계산
                        vec_OS = RotSatPos(SatPos, t_delay); % 지구자전효과 고려
                    end
                    %% DD 계산치 생성 파트 - 거리계산치를 각각 계산한 다음 DD 계산치 계산
                    vec_BsRS = vec_RS - TruePosBs_PP';  com_BsRS = norm(vec_BsRS);
                    vec_RvRS = vec_RS - x(1:3);    com_RvRS = norm(vec_RvRS);
                    vec_BsOS = vec_OS - TruePosBs_PP';  com_BsOS = norm(vec_BsOS);
                    vec_RvOS = vec_OS - x(1:3);    com_RvOS = norm(vec_RvOS);
                    com = (com_BsRS - com_RvRS) - (com_BsOS - com_RvOS);
                    y = obs - com;
                    Y(cnt2,1) = y;
                    
                    if Iter == 1
                        %% 각 위성 az, el 저장
                        [azBsRS,elBsRS] = xyz2azel(vec_BsRS, LatBs, LonBs);
                        azelBsRS(j,:) = [azBsRS, elBsRS, RS];
                        [azBsOS,elBsOS] = xyz2azel(vec_BsOS, LatBs, LonBs);
                        [azRvOS,elRvOS] = xyz2azel(vec_RvOS, LatRv, LonRv);
                        azelBsOS(cnt,:) = [azBsOS, elBsOS, OS];
                        %% Weighting
                        SNRRS(j,:) = [snr_Bs_RS, snr_Rv_RS,gs,RS];
                        SNROS(cnt,:) = [snr_Bs_OS, snr_Rv_OS,gs,OS];
                        DDel(cnt,:) = [elBsRS, elBsOS,gs,OS];
                        cnt = cnt + 1;
                        el = elBsOS;
                    end
                    if el >= eleCut
                        % 가중치 계산 파트
                        RS_snr = [snr_Bs_RS, snr_Rv_RS,gs,RS];
                        OS_snr = [snr_Bs_OS, snr_Rv_OS,gs,OS];
                        if cnt2 == 1
                            elsnr_Bs(cnt2,1:2) = [elBsRS, snr_Bs_RS];
                            elsnr_Rv(cnt2,1:2) = [elBsRS, snr_Rv_OS];
                        end
                        elsnr_Bs(cnt2+1,1:2) = [elBsOS, snr_Bs_OS];
                        elsnr_Rv(cnt2+1,1:2) = [elRvOS, snr_Rv_OS];
                        % H 행렬 계산 파트
                        H(cnt2,1) = vec_RvRS(1)/com_RvRS - vec_RvOS(1)/com_RvOS;
                        H(cnt2,2) = vec_RvRS(2)/com_RvRS - vec_RvOS(2)/com_RvOS;
                        H(cnt2,3) = vec_RvRS(3)/com_RvRS - vec_RvOS(3)/com_RvOS;
                        
                        cnt2= cnt2 + 1;
                        %                 HTH = HTH + H'*W*H;
                        %                 HTy = HTy + H'*W*y;
                    end
                    
                end
            elseif SYS == 4
                                %% 기준위성 좌표 먼저 계산 - Bs 기준
                icol = PickEPH(eph, RS ,gs);
                t_delay = QM_Bs_1e(find(QM_Bs_1e(:,2) == RS),4)/CCC;
                t_delay = (QM_Rv_1e(find(QM_Rv_1e(:,2) == RS),4)/(CCC+t_delay))/2;
                tc = gs - t_delay;
                vec_RS = GetSatPosNC_GC(eph, icol, tc);
                vec_RS = RotSatPos(vec_RS, t_delay);                                    % 지구 자전효과 고려
                %% GPS, BDS or GPS/BDS DD part
                for kS = 1:NoSats
                    if kS == indxRS || SatsEl(kS, 3) < eleCut
                        if SatsEl(kS, 3) < eleCut
                            NoSatsUsed = NoSatsUsed - 1;
                            if SatsEl(kS, 3) < 200
                                NoGPSSatsUsed = NoGPSSatsUsed - 1;
                            elseif SatsEl(kS, 3) > 200 & SatsEl(kS, 3) < 300
                                NoBDSSatsUsed = NoBDSSatsUsed - 1;
                            elseif SatsEl(kS, 3) > 300
                                NoGLOSatsUsed = NoGLOSatsUsed - 1;
                            end
                            disp([RS SatsEl(kS, 2)])
                        end
                        continue
                    end
                    %% DD 관측치 생성 파트
                    OS = Sats(kS,1);
                    OtherSats(kS,1) = OS;
                    snr_Bs_OS = QM_Bs_snr_1e(find(QM_Bs_snr_1e(:,2) == OS),4);           % 베이스, GPS 기준 위성 snr 값
                    snr_Rv_OS = QM_Rv_snr_1e(find(QM_Rv_snr_1e(:,2) == OS),4);           % 베이스, GPS 기준 위성 snr 값
                    
                    L1_Bs_OS = QM_Bs_L1_1e(find(QM_Bs_L1_1e(:,2) == OS),4);           % 베이스, GPS 기준 위성 snr 값
                    L1_Rv_OS = QM_Rv_L1_1e(find(QM_Rv_L1_1e(:,2) == OS),4);           % 베이스, GPS 기준 위성 snr 값
                    
                    obs_BsRS = QM_Bs_1e(find(QM_Bs_1e(:, 2) == RS), 4);
                    obs_RvRS = QM_Rv_1e(find(QM_Rv_1e(:, 2) == RS), 4);
                    obs_BsOS = QM_Bs_1e(find(QM_Bs_1e(:, 2) == OS), 4);
                    obs_RvOS = QM_Rv_1e(find(QM_Rv_1e(:, 2) == OS), 4);
                    
                    obs = (obs_BsRS - obs_RvRS) - (obs_BsOS - obs_RvOS);
                    obs2 = (L1_Bs_RS - L1_Rv_RS) - (L1_Bs_OS - L1_Rv_OS);
                    
                    OBS(cnt,1:2) = [(obs_BsRS - obs_RvRS) (obs_BsOS - obs_RvOS)];
                    %                         OBS2(cnt,1:2) = [(L1_Bs_RS - L1_Rv_RS) (L1_Bs_OS - L1_Rv_OS)];
                    %% DD 계산치 생성 파트 - 기타위성 좌표 계산(기준위성 좌표는 이미 계산 완료)
                    icol = PickEPH(eph, OS, gs);
                    t_delay = QM_Bs_1e(find(QM_Bs_1e(:,2) == OS),4)/CCC;
                    t_delay = (QM_Rv_1e(find(QM_Rv_1e(:,2) == OS),4)/(CCC+t_delay))/2;
                    tc = gs - t_delay;
                    vec_OS = GetSatPosNC_GC(eph, icol, tc);
                    vec_OS = RotSatPos(vec_OS, t_delay);
                    %% DD 계산치 생성 파트 - 거리계산치를 각각 계산한 다음 DD 계산치 계산
                    vec_BsRS = vec_RS - TruePosBs_PP';  com_BsRS = norm(vec_BsRS);
                    vec_RvRS = vec_RS - x(1:3);    com_RvRS = norm(vec_RvRS);
                    vec_BsOS = vec_OS - TruePosBs_PP';  com_BsOS = norm(vec_BsOS);
                    vec_RvOS = vec_OS - x(1:3);    com_RvOS = norm(vec_RvOS);
                    %                     vec_BsRS = TruePosBs_PP' - vec_RS;  com_BsRS = norm(vec_BsRS);
                    %                     vec_RvRS = x(1:3) - vec_RS;    com_RvRS = norm(vec_RvRS);
                    %                     vec_BsOS = TruePosBs_PP' - vec_OS;  com_BsOS = norm(vec_BsOS);
                    %                     vec_RvOS = x(1:3) - vec_OS;    com_RvOS = norm(vec_RvOS);
                    com = (com_BsRS - com_RvRS) - (com_BsOS - com_RvOS);
                    y = obs - com;
                    Y(cnt2,1) = y;
                    
                    if Iter == 1
                        %% 각 위성 az, el 저장
                        [azBsRS,elBsRS] = xyz2azel(vec_BsRS, LatBs, LonBs);
                        azelBsRS(j,:) = [azBsRS, elBsRS, RS];
                        [azBsOS,elBsOS] = xyz2azel(vec_BsOS, LatBs, LonBs);
                        [azRvOS,elRvOS] = xyz2azel(vec_RvOS, LatRv, LonRv);
                        azelBsOS(cnt,:) = [azBsOS, elBsOS, OS];
                        %% Weighting
                        SNRRS(j,:) = [snr_Bs_RS, snr_Rv_RS,gs,RS];
                        SNROS(cnt,:) = [snr_Bs_OS, snr_Rv_OS,gs,OS];
                        DDel(cnt,:) = [elBsRS, elBsOS,gs,OS];
                        cnt = cnt + 1;
                        el = elBsOS;
                    end
                    if el >= eleCut
                        % 가중치 계산 파트
                        RS_snr = [snr_Bs_RS, snr_Rv_RS,gs,RS];
                        OS_snr = [snr_Bs_OS, snr_Rv_OS,gs,OS];
                        if cnt2 == 1
                            elsnr_Bs(cnt2,1:2) = [elBsRS, snr_Bs_RS];
                            elsnr_Rv(cnt2,1:2) = [elBsRS, snr_Rv_OS];
                        end
                        elsnr_Bs(cnt2+1,1:2) = [elBsOS, snr_Bs_OS];
                        elsnr_Rv(cnt2+1,1:2) = [elRvOS, snr_Rv_OS];
                        H(cnt2,1) = vec_RvRS(1)/com_RvRS - vec_RvOS(1)/com_RvOS;
                        H(cnt2,2) = vec_RvRS(2)/com_RvRS - vec_RvOS(2)/com_RvOS;
                        H(cnt2,3) = vec_RvRS(3)/com_RvRS - vec_RvOS(3)/com_RvOS;
                        if OS < 200
                            H(cnt2,4) = 0;
                        else
                            H(cnt2,4) = 1;
                        end
                        cnt2= cnt2 + 1;
                    end
                    
                end
            elseif SYS == 5 | SYS == 6
                %% 기준위성 좌표 먼저 계산 - Bs 기준
                if RS < 300
                    icol = PickEPH(eph, RS ,gs);
                    t_delay = QM_Bs_1e(find(QM_Bs_1e(:,2) == RS),4)/CCC;
                    t_delay = (QM_Rv_1e(find(QM_Rv_1e(:,2) == RS),4)/(CCC+t_delay))/2;
                    %                     tc = gs - STT;
                    tc = gs - t_delay;
                    vec_RS = GetSatPosNC_GC(eph, icol, tc);
                    vec_RS = RotSatPos(vec_RS, t_delay);                                    % 지구 자전효과 고려
                elseif RS > 300
                    icol = PickEPH_GLO2(ephGLO, RS-300, gs);
                    t_delay = QM_Bs_1e(find(QM_Bs_1e(:,2) == RS),4)/CCC;
                    t_delay = (QM_Rv_1e(find(QM_Rv_1e(:,2) == RS),4)/(CCC+t_delay))/2;
                    tc = gs - t_delay;
                    tc = tc - LeapSec - TauC;
                    [vec_RS,SatVel,SatLS] = GetSatPosGLO(ephGLO,icol,tc,deltat); % 방송궤도력으로 위성위치 계산
                    vec_RS = RotSatPos(vec_RS, t_delay); % 지구자전효과 고려
                end
                %% GPS, BDS or GPS/BDS DD part
                for kS = 1:NoSats
                    if kS == indxRS || SatsEl(kS, 3) < eleCut
                        if SatsEl(kS, 3) < eleCut
                            NoSatsUsed = NoSatsUsed - 1;
                            if SatsEl(kS, 3) < 200
                                NoGPSSatsUsed = NoGPSSatsUsed - 1;
                            elseif SatsEl(kS, 3) > 200 & SatsEl(kS, 3) < 300
                                NoBDSSatsUsed = NoBDSSatsUsed - 1;
                            elseif SatsEl(kS, 3) > 300
                                NoGLOSatsUsed = NoGLOSatsUsed - 1;
                            end
                            disp([RS SatsEl(kS, 2)])
                        end
                        continue
                    end
                    %% DD 관측치 생성 파트
                    OS = Sats(kS,1);
                    OtherSats(kS,1) = OS;
                    snr_Bs_OS = QM_Bs_snr_1e(find(QM_Bs_snr_1e(:,2) == OS),4);           % 베이스, GPS 기준 위성 snr 값
                    snr_Rv_OS = QM_Rv_snr_1e(find(QM_Rv_snr_1e(:,2) == OS),4);           % 베이스, GPS 기준 위성 snr 값
                    
                    L1_Bs_OS = QM_Bs_L1_1e(find(QM_Bs_L1_1e(:,2) == OS),4);           % 베이스, GPS 기준 위성 snr 값
                    L1_Rv_OS = QM_Rv_L1_1e(find(QM_Rv_L1_1e(:,2) == OS),4);           % 베이스, GPS 기준 위성 snr 값
                    
                    obs_BsRS = QM_Bs_1e(find(QM_Bs_1e(:, 2) == RS), 4);
                    obs_RvRS = QM_Rv_1e(find(QM_Rv_1e(:, 2) == RS), 4);
                    obs_BsOS = QM_Bs_1e(find(QM_Bs_1e(:, 2) == OS), 4);
                    obs_RvOS = QM_Rv_1e(find(QM_Rv_1e(:, 2) == OS), 4);
                    
                    obs = (obs_BsRS - obs_RvRS) - (obs_BsOS - obs_RvOS);
                    obs2 = (L1_Bs_RS - L1_Rv_RS) - (L1_Bs_OS - L1_Rv_OS);
                    
                    OBS(cnt,1:2) = [(obs_BsRS - obs_RvRS) (obs_BsOS - obs_RvOS)];
                    %% DD 계산치 생성 파트 - 기타위성 좌표 계산(기준위성 좌표는 이미 계산 완료)
                    if OS < 300
                        icol = PickEPH(eph, OS, gs);
                        t_delay = QM_Bs_1e(find(QM_Bs_1e(:,2) == OS),4)/CCC;
                        t_delay = (QM_Rv_1e(find(QM_Rv_1e(:,2) == OS),4)/(CCC+t_delay))/2;
                        tc = gs - t_delay;
                        vec_OS = GetSatPosNC_GC(eph, icol, tc);
                        vec_OS = RotSatPos(vec_OS, t_delay);
                        %% DD 계산치 생성 파트 - 거리계산치를 각각 계산한 다음 DD 계산치 계산
                        vec_BsRS = vec_RS - TruePosBs_PP';  com_BsRS = norm(vec_BsRS);
                        vec_RvRS = vec_RS - x(1:3);    com_RvRS = norm(vec_RvRS);
                        vec_BsOS = vec_OS - TruePosBs_PP';  com_BsOS = norm(vec_BsOS);
                        vec_RvOS = vec_OS - x(1:3);    com_RvOS = norm(vec_RvOS);
                        com = (com_BsRS - com_RvRS) - (com_BsOS - com_RvOS);
                        y = obs - com;
                        Y(cnt2,1) = y;
                        
                        if Iter == 1
                            %% 각 위성 az, el 저장
                            [azBsRS,elBsRS] = xyz2azel(vec_BsRS, LatBs, LonBs);
                            azelBsRS(j,:) = [azBsRS, elBsRS, RS];
                            [azBsOS,elBsOS] = xyz2azel(vec_BsOS, LatBs, LonBs);
                            [azRvOS,elRvOS] = xyz2azel(vec_RvOS, LatRv, LonRv);
                            azelBsOS(cnt,:) = [azBsOS, elBsOS, OS];
                            %% Weighting
                            SNRRS(j,:) = [snr_Bs_RS, snr_Rv_RS,gs,RS];
                            SNROS(cnt,:) = [snr_Bs_OS, snr_Rv_OS,gs,OS];
                            DDel(cnt,:) = [elBsRS, elBsOS,gs,OS];
                            cnt = cnt + 1;
                            el = elBsOS;
                        end
                        if el >= eleCut
                            % 가중치 계산 파트
                            RS_snr = [snr_Bs_RS, snr_Rv_RS,gs,RS];
                            OS_snr = [snr_Bs_OS, snr_Rv_OS,gs,OS];
                            if cnt2 == 1
                                elsnr_Bs(cnt2,1:2) = [elBsRS, snr_Bs_RS];
                                elsnr_Rv(cnt2,1:2) = [elBsRS, snr_Rv_OS];
                            end
                            elsnr_Bs(cnt2+1,1:2) = [elBsOS, snr_Bs_OS];
                            elsnr_Rv(cnt2+1,1:2) = [elRvOS, snr_Rv_OS];
                            
                            H(cnt2,1) = vec_RvRS(1)/com_RvRS - vec_RvOS(1)/com_RvOS;
                            H(cnt2,2) = vec_RvRS(2)/com_RvRS - vec_RvOS(2)/com_RvOS;
                            H(cnt2,3) = vec_RvRS(3)/com_RvRS - vec_RvOS(3)/com_RvOS;
                            if SYS == 4
                                if OS < 200
                                    H(cnt2,4) = 0;
                                else
                                    H(cnt2,4) = 1;
                                end
                            else
                                if OS < 300
                                    H(cnt2,4) = 0;
                                else
                                    H(cnt2,4) = 1;
                                end
                            end
                            cnt2= cnt2 + 1;
                        end
                    elseif OS > 300
                        icol = PickEPH_GLO2(ephGLO, OS-300, gs);
                        icolArr(kS,1)=OS; icolArr(kS,2)=icol;
                        TauN = ephGLO(icol,12); % : tau & gamma - 시계오차 보정에 필요
                        GammaN = ephGLO(icol,13);
                        ch_num = ephGLO(icol,16); % : channel number - 전리층 보정에 필요
                        jcol = find(SatPosArr_before(:,1)==OS);
                        if isempty(jcol); jcol=0;end
                        jcol=jcol(1);
                        %% 위성위치 계산 파트
                        if icol == icolArr_before(icolArr_before(:,1)==OS,2) % icol 변화 없으면 이전 epoch으로 계산
                            %                 disp(1)
                            t_delay = QM_Bs_1e(find(QM_Bs_1e(:,2) == OS),4)/CCC;
                            t_delay = (QM_Rv_1e(find(QM_Rv_1e(:,2) == OS),4)/(CCC+t_delay))/2;
                            tc = gs - t_delay;
                            tc = tc - LeapSec - TauC;
                            [SatPos,SatVel,SatLS] = GetSatPosGLO_ver3(SatPosArr_before,jcol,tc,deltat); % 이전 epoch값으로 위성위치 계산
                            vec_OS = RotSatPos(SatPos, t_delay); % 지구자전효과 고려
                        else % icol 변화시(새위성 또는 새로운 방송궤도력) 방송궤도력으로 계산
                            %                 disp(2)
                            t_delay = QM_Bs_1e(find(QM_Bs_1e(:,2) == OS),4)/CCC;
                            t_delay = (QM_Rv_1e(find(QM_Rv_1e(:,2) == OS),4)/(CCC+t_delay))/2;
                            tc = gs - t_delay;
                            tc = tc - LeapSec - TauC;
                            [SatPos,SatVel,SatLS] = GetSatPosGLO(ephGLO,icol,tc,deltat); % 방송궤도력으로 위성위치 계산
                            vec_OS = RotSatPos(SatPos, t_delay); % 지구자전효과 고려
                        end
                        %% DD 계산치 생성 파트 - 거리계산치를 각각 계산한 다음 DD 계산치 계산
                        vec_BsRS = vec_RS - TruePosBs_PP';  com_BsRS = norm(vec_BsRS);
                        vec_RvRS = vec_RS - x(1:3);    com_RvRS = norm(vec_RvRS);
                        vec_BsOS = vec_OS - TruePosBs_PP';  com_BsOS = norm(vec_BsOS);
                        vec_RvOS = vec_OS - x(1:3);    com_RvOS = norm(vec_RvOS);
                        com = (com_BsRS - com_RvRS) - (com_BsOS - com_RvOS);
                        y = obs - com;
                        Y(cnt2,1) = y;
                        
                        if Iter == 1
                            %% 각 위성 az, el 저장
                            [azBsRS,elBsRS] = xyz2azel(vec_BsRS, LatBs, LonBs);
                            azelBsRS(j,:) = [azBsRS, elBsRS, RS];
                            [azBsOS,elBsOS] = xyz2azel(vec_BsOS, LatBs, LonBs);
                            [azRvOS,elRvOS] = xyz2azel(vec_RvOS, LatRv, LonRv);
                            azelBsOS(cnt,:) = [azBsOS, elBsOS, OS];
                            %% Weighting
                            SNRRS(j,:) = [snr_Bs_RS, snr_Rv_RS,gs,RS];
                            SNROS(cnt,:) = [snr_Bs_OS, snr_Rv_OS,gs,OS];
                            DDel(cnt,:) = [elBsRS, elBsOS,gs,OS];
                            cnt = cnt + 1;
                            el = elBsOS;
                        end
                        if el >= eleCut
                            % 가중치 계산 파트
                            RS_snr = [snr_Bs_RS, snr_Rv_RS,gs,RS];
                            OS_snr = [snr_Bs_OS, snr_Rv_OS,gs,OS];
                            w = 1;
                            if cnt2 == 1
                                elsnr_Bs(cnt2,1:2) = [elBsRS, snr_Bs_RS];
                                elsnr_Rv(cnt2,1:2) = [elBsRS, snr_Rv_OS];
                            end
                            elsnr_Bs(cnt2+1,1:2) = [elBsOS, snr_Bs_OS];
                            elsnr_Rv(cnt2+1,1:2) = [elRvOS, snr_Rv_OS];
                            % H 행렬 계산 파트
                            H(cnt2,1) = vec_RvRS(1)/com_RvRS - vec_RvOS(1)/com_RvOS;
                            H(cnt2,2) = vec_RvRS(2)/com_RvRS - vec_RvOS(2)/com_RvOS;
                            H(cnt2,3) = vec_RvRS(3)/com_RvRS - vec_RvOS(3)/com_RvOS;
                            H(cnt2,4) = 1;
                            
                            
                            cnt2= cnt2 + 1;
                        end
                        
                    end
                end
            elseif SYS == 7
                %% 기준위성 좌표 먼저 계산 - Bs 기준
                if RS < 300
                    icol = PickEPH(eph, RS ,gs);
                    t_delay = QM_Bs_1e(find(QM_Bs_1e(:,2) == RS),4)/CCC;
                    t_delay = (QM_Rv_1e(find(QM_Rv_1e(:,2) == RS),4)/(CCC+t_delay))/2;
                    %                     tc = gs - STT;
                    tc = gs - t_delay;
                    vec_RS = GetSatPosNC_GC(eph, icol, tc);
                    vec_RS = RotSatPos(vec_RS, t_delay);                                    % 지구 자전효과 고려
                elseif RS > 300
                    icol = PickEPH_GLO2(ephGLO, RS-300, gs);
                    t_delay = QM_Bs_1e(find(QM_Bs_1e(:,2) == RS),4)/CCC;
                    t_delay = (QM_Rv_1e(find(QM_Rv_1e(:,2) == RS),4)/(CCC+t_delay))/2;
                    tc = gs - t_delay;
                    tc = tc - LeapSec - TauC;
                    [vec_RS,SatVel,SatLS] = GetSatPosGLO(ephGLO,icol,tc,deltat); % 방송궤도력으로 위성위치 계산
                    vec_RS = RotSatPos(vec_RS, t_delay); % 지구자전효과 고려
                end
                %% GPS, BDS or GPS/BDS DD part
                for kS = 1:NoSats
                    if kS == indxRS || SatsEl(kS, 3) < eleCut
                        if SatsEl(kS, 3) < eleCut
                            NoSatsUsed = NoSatsUsed - 1;
                            if SatsEl(kS, 3) < 200
                                NoGPSSatsUsed = NoGPSSatsUsed - 1;
                            elseif SatsEl(kS, 3) > 200 & SatsEl(kS, 3) < 300
                                NoBDSSatsUsed = NoBDSSatsUsed - 1;
                            elseif SatsEl(kS, 3) > 300
                                NoGLOSatsUsed = NoGLOSatsUsed - 1;
                            end
                            disp([RS SatsEl(kS, 2)])
                        end
                        continue
                    end
                    %% DD 관측치 생성 파트
                    OS = Sats(kS,1);
                    OtherSats(kS,1) = OS;
                    snr_Bs_OS = QM_Bs_snr_1e(find(QM_Bs_snr_1e(:,2) == OS),4);           % 베이스, GPS 기준 위성 snr 값
                    snr_Rv_OS = QM_Rv_snr_1e(find(QM_Rv_snr_1e(:,2) == OS),4);           % 베이스, GPS 기준 위성 snr 값
                    
                    L1_Bs_OS = QM_Bs_L1_1e(find(QM_Bs_L1_1e(:,2) == OS),4);           % 베이스, GPS 기준 위성 snr 값
                    L1_Rv_OS = QM_Rv_L1_1e(find(QM_Rv_L1_1e(:,2) == OS),4);           % 베이스, GPS 기준 위성 snr 값
                    
                    obs_BsRS = QM_Bs_1e(find(QM_Bs_1e(:, 2) == RS), 4);
                    obs_RvRS = QM_Rv_1e(find(QM_Rv_1e(:, 2) == RS), 4);
                    obs_BsOS = QM_Bs_1e(find(QM_Bs_1e(:, 2) == OS), 4);
                    obs_RvOS = QM_Rv_1e(find(QM_Rv_1e(:, 2) == OS), 4);
                    
                    obs = (obs_BsRS - obs_RvRS) - (obs_BsOS - obs_RvOS);
                    obs2 = (L1_Bs_RS - L1_Rv_RS) - (L1_Bs_OS - L1_Rv_OS);
                    
                    OBS(cnt,1:2) = [(obs_BsRS - obs_RvRS) (obs_BsOS - obs_RvOS)];
                    %% DD 계산치 생성 파트 - 기타위성 좌표 계산(기준위성 좌표는 이미 계산 완료)
                    if OS < 300
                        icol = PickEPH(eph, OS, gs);
                        t_delay = QM_Bs_1e(find(QM_Bs_1e(:,2) == OS),4)/CCC;
                        t_delay = (QM_Rv_1e(find(QM_Rv_1e(:,2) == OS),4)/(CCC+t_delay))/2;
                        tc = gs - t_delay;
                        vec_OS = GetSatPosNC_GC(eph, icol, tc);
                        vec_OS = RotSatPos(vec_OS, t_delay);
                        %% DD 계산치 생성 파트 - 거리계산치를 각각 계산한 다음 DD 계산치 계산
                        vec_BsRS = vec_RS - TruePosBs_PP';  com_BsRS = norm(vec_BsRS);
                        vec_RvRS = vec_RS - x(1:3);    com_RvRS = norm(vec_RvRS);
                        vec_BsOS = vec_OS - TruePosBs_PP';  com_BsOS = norm(vec_BsOS);
                        vec_RvOS = vec_OS - x(1:3);    com_RvOS = norm(vec_RvOS);
                        com = (com_BsRS - com_RvRS) - (com_BsOS - com_RvOS);
                        y = obs - com;
                        Y(cnt2,1) = y;
                        
                        if Iter == 1
                            %% 각 위성 az, el 저장
                            [azBsRS,elBsRS] = xyz2azel(vec_BsRS, LatBs, LonBs);
                            azelBsRS(j,:) = [azBsRS, elBsRS, RS];
                            [azBsOS,elBsOS] = xyz2azel(vec_BsOS, LatBs, LonBs);
                            [azRvOS,elRvOS] = xyz2azel(vec_RvOS, LatRv, LonRv);
                            azelBsOS(cnt,:) = [azBsOS, elBsOS, OS];
                            %% Weighting
                            SNRRS(j,:) = [snr_Bs_RS, snr_Rv_RS,gs,RS];
                            SNROS(cnt,:) = [snr_Bs_OS, snr_Rv_OS,gs,OS];
                            DDel(cnt,:) = [elBsRS, elBsOS,gs,OS];
                            cnt = cnt + 1;
                            el = elBsOS;
                        end
                        if el >= eleCut
                            % 가중치 계산 파트
                            RS_snr = [snr_Bs_RS, snr_Rv_RS,gs,RS];
                            OS_snr = [snr_Bs_OS, snr_Rv_OS,gs,OS];
                            if cnt2 == 1
                                elsnr_Bs(cnt2,1:2) = [elBsRS, snr_Bs_RS];
                                elsnr_Rv(cnt2,1:2) = [elBsRS, snr_Rv_OS];
                            end
                            elsnr_Bs(cnt2+1,1:2) = [elBsOS, snr_Bs_OS];
                            elsnr_Rv(cnt2+1,1:2) = [elRvOS, snr_Rv_OS];
                            
                            H(cnt2,1) = vec_RvRS(1)/com_RvRS - vec_RvOS(1)/com_RvOS;
                            H(cnt2,2) = vec_RvRS(2)/com_RvRS - vec_RvOS(2)/com_RvOS;
                            H(cnt2,3) = vec_RvRS(3)/com_RvRS - vec_RvOS(3)/com_RvOS;
                            if OS < 200
                                H(cnt2,4:5) = [0 0];
                            elseif OS > 200 & OS < 300
                                H(cnt2,4:5) = [1 0];
                            end
                            
                            cnt2= cnt2 + 1;
                        end
                    elseif OS > 300
                        icol = PickEPH_GLO2(ephGLO, OS-300, gs);
                        icolArr(kS,1)=OS; icolArr(kS,2)=icol;
                        TauN = ephGLO(icol,12); % : tau & gamma - 시계오차 보정에 필요
                        GammaN = ephGLO(icol,13);
                        ch_num = ephGLO(icol,16); % : channel number - 전리층 보정에 필요
                        jcol = find(SatPosArr_before(:,1)==OS);
                        if isempty(jcol); jcol=0;end
                        jcol=jcol(1);
                        %% 위성위치 계산 파트
                        if icol == icolArr_before(icolArr_before(:,1)==OS,2) % icol 변화 없으면 이전 epoch으로 계산
                            %                 disp(1)
                            t_delay = QM_Bs_1e(find(QM_Bs_1e(:,2) == OS),4)/CCC;
                            t_delay = (QM_Rv_1e(find(QM_Rv_1e(:,2) == OS),4)/(CCC+t_delay))/2;
                            tc = gs - t_delay;
                            tc = tc - LeapSec - TauC;
                            [SatPos,SatVel,SatLS] = GetSatPosGLO_ver3(SatPosArr_before,jcol,tc,deltat); % 이전 epoch값으로 위성위치 계산
                            vec_OS = RotSatPos(SatPos, t_delay); % 지구자전효과 고려
                        else % icol 변화시(새위성 또는 새로운 방송궤도력) 방송궤도력으로 계산
                            %                 disp(2)
                            t_delay = QM_Bs_1e(find(QM_Bs_1e(:,2) == OS),4)/CCC;
                            t_delay = (QM_Rv_1e(find(QM_Rv_1e(:,2) == OS),4)/(CCC+t_delay))/2;
                            tc = gs - t_delay;
                            tc = tc - LeapSec - TauC;
                            [SatPos,SatVel,SatLS] = GetSatPosGLO(ephGLO,icol,tc,deltat); % 방송궤도력으로 위성위치 계산
                            vec_OS = RotSatPos(SatPos, t_delay); % 지구자전효과 고려
                        end
                        %% DD 계산치 생성 파트 - 거리계산치를 각각 계산한 다음 DD 계산치 계산
                        vec_BsRS = vec_RS - TruePosBs_PP';  com_BsRS = norm(vec_BsRS);
                        vec_RvRS = vec_RS - x(1:3);    com_RvRS = norm(vec_RvRS);
                        vec_BsOS = vec_OS - TruePosBs_PP';  com_BsOS = norm(vec_BsOS);
                        vec_RvOS = vec_OS - x(1:3);    com_RvOS = norm(vec_RvOS);
                        com = (com_BsRS - com_RvRS) - (com_BsOS - com_RvOS);
                        y = obs - com;
                        Y(cnt2,1) = y;
                        
                        if Iter == 1
                            %% 각 위성 az, el 저장
                            [azBsRS,elBsRS] = xyz2azel(vec_BsRS, LatBs, LonBs);
                            azelBsRS(j,:) = [azBsRS, elBsRS, RS];
                            [azBsOS,elBsOS] = xyz2azel(vec_BsOS, LatBs, LonBs);
                            [azRvOS,elRvOS] = xyz2azel(vec_RvOS, LatRv, LonRv);
                            azelBsOS(cnt,:) = [azBsOS, elBsOS, OS];
                            %% Weighting
                            SNRRS(j,:) = [snr_Bs_RS, snr_Rv_RS,gs,RS];
                            SNROS(cnt,:) = [snr_Bs_OS, snr_Rv_OS,gs,OS];
                            DDel(cnt,:) = [elBsRS, elBsOS,gs,OS];
                            cnt = cnt + 1;
                            el = elBsOS;
                        end
                        if el >= eleCut
                            % 가중치 계산 파트
                            RS_snr = [snr_Bs_RS, snr_Rv_RS,gs,RS];
                            OS_snr = [snr_Bs_OS, snr_Rv_OS,gs,OS];
                            w = 1;
                            if cnt2 == 1
                                elsnr_Bs(cnt2,1:2) = [elBsRS, snr_Bs_RS];
                                elsnr_Rv(cnt2,1:2) = [elBsRS, snr_Rv_OS];
                            end
                            elsnr_Bs(cnt2+1,1:2) = [elBsOS, snr_Bs_OS];
                            elsnr_Rv(cnt2+1,1:2) = [elRvOS, snr_Rv_OS];
                            % H 행렬 계산 파트
                            H(cnt2,1) = vec_RvRS(1)/com_RvRS - vec_RvOS(1)/com_RvOS;
                            H(cnt2,2) = vec_RvRS(2)/com_RvRS - vec_RvOS(2)/com_RvOS;
                            H(cnt2,3) = vec_RvRS(3)/com_RvRS - vec_RvOS(3)/com_RvOS;
                            H(cnt2,4:5) = [0 1];
                            
                            
                            cnt2= cnt2 + 1;
                        end
                        
                    end
                end
                
            end
            
            W = DD_cofactor_matrix(elsnr_Bs(:,1), elsnr_Rv(:,1), elsnr_Bs(:,2), elsnr_Rv(:,2), 4);
            HTH = H(1:cnt2-1,:)'*inv(W(1:cnt2-1,1:cnt2-1))*H(1:cnt2-1,:);
            HTy = H(1:cnt2-1,:)'*inv(W(1:cnt2-1,1:cnt2-1))*Y(1:cnt2-1,:);
            
            xhat = inv(HTH) * HTy;
            x = x + xhat;
            XHAT(j,1:2) =[gs, norm(xhat)];
            if SYS < 4
                if norm(xhat) < EpsStop;
                    nEst = nEst + 1;
                    estm(nEst,1) =gs;
                    estm(nEst,2:4) =x;
                    estm(nEst,5) = length(GPS(:,1));        % GPS 관측 위성 수
                    estm(nEst,6) = length(BDS(:,1));        % BDS 관측 위성 수
                    estm(nEst,7) = NoGPSSatsUsed;           % : snr issue 가 없을때
                    estm(nEst,8) = NoBDSSatsUsed;           % : snr issue 가 없을때
                    estm(nEst,9) = NoGLOSatsUsed;           % : snr issue 가 없을때
                    estm(nEst,10) = NoSatsUsed;           % : snr issue 가 없을때
                    toc;
                    break;
                end
            else
                if Iter == MaxIter
                    nEst = nEst + 1;
                    estm(nEst,1) =gs;
                    estm(nEst,2:4) =x(1:3);
                    estm(nEst,5) = length(GPS(:,1));        % GPS 관측 위성 수
                    estm(nEst,6) = length(BDS(:,1));        % BDS 관측 위성 수
                    estm(nEst,7) = NoGPSSatsUsed;           % : snr issue 가 없을때
                    estm(nEst,8) = NoBDSSatsUsed;           % : snr issue 가 없을때
                    estm(nEst,9) = NoGLOSatsUsed;           % : snr issue 가 없을때
                    estm(nEst,10) = NoSatsUsed;           % : snr issue 가 없을때
                end
            end
        end
    end
    
    %% rover's longi, lati
    rover_gd(nEst,1) = gs;
    rover_gd(nEst,2:4) = xyz2gd(estm(nEst,2:4)); % rover의 xyz값 gd 로 변환
    AppLat = rover_gd(nEst,2); AppLon = rover_gd(nEst,3);
    
end
estm_x = smooth(estm(:,2),0.1,'rloess');
estm_y = smooth(estm(:,3),0.1,'rloess');
estm_z = smooth(estm(:,4),0.1,'rloess');
estm_new = estm;
estm_new(:,2:4) = [estm_x, estm_y, estm_z];

azel_sum = [azelBsOS; azelBsRS];
snr_sum = [SNROS;SNRRS];
data = [azel_sum(:,3), azel_sum(:,2), azel_sum(:,1), snr_sum(:,1)];
%% 측위오차 분석 & 그래프 작성
% estm = estm(650:end, :);
% Base = Base(650:end, :);
if Dynamic == 0
    [DDdXYZ, DDdNEV, DDdis, DDrms, DDstd] = PostErrorsDDstatic(estm, Base, Truedis,8, 12, SYS);    % 임의의 장소에서 이동 측위시
elseif Dynamic == 1
    [DDdXYZ, DDdXYZ_vrs, DDdNEV, DDdNEV_vrs, DDdis, DDrms, result] = PostErrorsDDkine(estm, Base, Base_vrs, Rover_vrs, -2, 2, SYS);    % 임의의 장소에서 이동 측위시
end

%% 각 위성별 SNR Plot
% DDPlotQM(QM11, QM22, 141, 'Base', 'Rover')

% skyplot = SkyPlot();
% for i=1:length(data)
% prn = data(i,1);
% az = data(i,2);
% el = data(i,3);
% SNR =data(i,4);
% skyplot.add(prn,az,el,SNR);
% end
% skyplot.PLOT();

if DOY == 25 & YY == 17
    [dXYZ, dNEV, result] = PostErrorskine(estm, [-3041241.74100000,4053944.14300000,3859873.64000000], SYS); 
else
    [dXYZ, dNEV, result] = PostErrorskine(estm, Rover_vrs, SYS);
end

