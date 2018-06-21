%% GPS & BDS 통합 측위 코드 (DGNSS)
% tic;
clear all; close all;
warning off;
%% 상수 정의
CCC = 299792458.;   % CCC = Speed of Light [m/s]
%% QM 파일, 날짜 핸들링
%-----------170116 테스트---------%
RTCM = load('PPS1_170116.t41');
DOY = 016;YY  = 17;
%% Navigation file load
FileNav = strcat('brdm', num2str(DOY,'%03d'), '0.', num2str(YY,'%02d'), 'p');   %: Navigation RINEX file
[eph, trashPRN, trashT]=ReadEPH_all(FileNav);
% load('PTH_PP_dop.mat');

for LR=7:8
    % for LR=1:24
    LR
    tic;
    switch LR
        %----세단 Type A----%
        case 1
            clear -global PAST_QM
            FileQM='test1_L';
        case 2
            FileQM='test1_R';
            clear -global PAST_QM
            %----세단 Type B----%
        case 3
            clear -global PAST_QM
            FileQM='test2_1_L';
        case 4
            clear -global PAST_QM
            FileQM='test2_1_R';
        case 5
            clear -global PAST_QM
            FileQM='test2_2_L';
        case 6
            clear -global PAST_QM
            FileQM='test2_2_R';
            %----SUV Type A----%
        case 7
            clear -global PAST_QM
            FileQM='test3_L';
        case 8
            clear -global PAST_QM
            FileQM='test3_R';
            %----SUV Type B----%
        case 9
            clear -global PAST_QM
            FileQM='test4_1_L';
        case 10
            clear -global PAST_QM
            FileQM='test4_1_R';
        case 11
            clear -global PAST_QM
            FileQM='test4_2_L';
        case 12
            clear -global PAST_QM
            FileQM='test4_2_R';
        case 13
            clear -global PAST_QM
            FileQM='test4_3_L';
        case 14
            clear -global PAST_QM
            FileQM='test4_3_R';
        case 15
            FileQM='testVT_L';
        case 16
            FileQM='testVT_R';
    end
    %-----------------------------------------------------%
    [arrQM, FinalPRNs, FinalTTs] = ReadQM(FileQM);
    arrQM = arrQM(find(arrQM(:,3) < 300),:);
    
    %% 이상 위성 제거
    switch LR
        case 1
            arrQM = arrQM(find(arrQM(:,1) > 1000),:);
            idx_wrong = [min(find(arrQM(:,1) == 105219)),max(find(arrQM(:,1) == 105343))];      % test1_L
            arrQM_1st = arrQM(1:idx_wrong(1)-1,:);
            arrQM_2nd = arrQM(idx_wrong(2)+1:end,:);
            arrQM_10 = arrQM(idx_wrong(1):idx_wrong(2),:);
            arrQM_10 = arrQM_10(find(arrQM_10(:,2) ~= 20),:);       % test1_L
            arrQM = [arrQM_1st; arrQM_10; arrQM_2nd];
        case 4
            idx_wrong = [min(find(arrQM(:,1) == 106284)),max(find(arrQM(:,1) == 106430))];      % test2_1_R
            arrQM_1st = arrQM(1:idx_wrong(1)-1,:);
            arrQM_2nd = arrQM(idx_wrong(2)+1:end,:);
            arrQM_10 = arrQM(idx_wrong(1):idx_wrong(2),:);
            arrQM_10 = arrQM_10(find(arrQM_10(:,2) ~= 10),:);       % test2_1_R
            arrQM = [arrQM_1st; arrQM_10; arrQM_2nd];
        case 6
            idx_wrong = [min(find(arrQM(:,1) == 107650)),max(find(arrQM(:,1) == 107676))];      % test2_2_R
            arrQM_1st = arrQM(1:idx_wrong(1)-1,:);
            arrQM_2nd = arrQM(idx_wrong(2)+1:end,:);
            arrQM_10 = arrQM(idx_wrong(1):idx_wrong(2),:);
            arrQM_10 = arrQM_10(find(arrQM_10(:,2) ~= 32),:);       % test2_2_R
            arrQM = [arrQM_1st; arrQM_10; arrQM_2nd];
        case 9
            idx_wrong = [min(find(arrQM(:,1) == 110332)),max(find(arrQM(:,1) == 110362))];      % test4_1_L
            arrQM_1st = arrQM(1:idx_wrong(1)-1,:);
            arrQM_2nd = arrQM(idx_wrong(2)+1:end,:);
            arrQM_10 = arrQM(idx_wrong(1):idx_wrong(2),:);
            arrQM_10 = arrQM_10(find(arrQM_10(:,2) ~= 25),:);       % test4_1_R
            arrQM = [arrQM_1st; arrQM_10; arrQM_2nd];
        case 10
            idx_wrong = [min(find(arrQM(:,1) == 110054)),max(find(arrQM(:,1) == 110074))];      % test4_1_R
            arrQM_1st = arrQM(1:idx_wrong(1)-1,:);
            arrQM_2nd = arrQM(idx_wrong(2)+1:end,:);
            arrQM_10 = arrQM(idx_wrong(1):idx_wrong(2),:);
            arrQM_10 = arrQM_10(find(arrQM_10(:,2) ~= 25),:);       % test4_1_R
            arrQM = [arrQM_1st; arrQM_10; arrQM_2nd];
    end
    % break
    yyS = num2str(YY,'%02d');
    doyS = num2str(DOY,'%03d');
    [gw, gd] = ydoy2gwgd(YY, DOY); %: GPS WEEK 결정
    [bw, bd] = ydoy2bwbd(YY, DOY); %: BDS WEEK 결정
    
    %% 좌표 참값
    TruePos = [-3058799.61420451,4083265.35912516,3814946.87192938]; %170116 평택항 초기좌표
    
    %% 윤초 핸들링
    LeapSecBDS = 14;
    
    %% brdm파일에는 al/be가 없기 때문에 brdc파일에서 가져옴
    al=[0 0 0 0];   be=[0 0 0 0];
    
    %% 사용할 관측치 선정
    g_ObsType = 120; % gps C1
    g_ObsType_snr = 141;
    c_ObsType = 220; % bds C1
    c_ObsType_snr = 241;
    
    %% QM 준비
    QM = SelectQM_gc(arrQM, g_ObsType, c_ObsType);
    QM_snr = SelectQM_gc(arrQM, g_ObsType_snr, c_ObsType_snr);
    FinalPRNs = unique(QM(:,2));
    FinalTTs = unique(QM(:,1));
    FinalTTs = FinalTTs(find(FinalTTs(:,1) > 1000),:);
    % break
    
    %% Doppler smoothing을 위한 QM 핸들링
    QMd =qmHandle(arrQM);
    
    %% 좌표 초기치 설정 및 위경도 변환
    AppPos = TruePos;
    gd = xyz2gd(AppPos); AppLat = gd(1); AppLon = gd(2);
    
    %% 추정을 위한 매개변수 설정
    Maxiter = 10;
    EpsStop = 1e-4;
    ctr = 1; ctr2 = 1;
    eleCut = 15;
    x = [AppPos ctr ctr2]';
    %%
    NoEpochs = length(FinalTTs);
    estm = zeros(NoEpochs, 6);
    Scount = zeros(NoEpochs, 1);
    nEst = 0;
    for j = 1:NoEpochs % 60:65%
        gs = FinalTTs(j);
        x = [AppPos ctr ctr2]';
%         clear -global PAST_QM
        QMe_ = qmHandle(QMd.pickQM(gs,':',':'));
        QMe = dopplerSM(QMe_.getQM(),100);
        indexQM = find(QM(:,1) == gs);
        QM1e = QM(indexQM,:);
        QM1e(:,4) = QMe(find(QMe(:,3) < 320),4);
        QM1e_snr = QM_snr(indexQM,:);
        NoSats = length(QM1e);
        NoSatsUsed = 0;
        for iter = 1:Maxiter
            
            HTH = zeros(5,5);
            HTy = zeros(5,1);
                                   
            if NoSats <= 6
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
                STT = GetSTTbrdm(gs, eph, icol, x(1:3)); % 신호전달시간 계산
                if prn > 100 && prn < 200
                    tc = gs - STT ;             % 신호전달시간 보정
                elseif prn > 200 && prn < 300
                    tc = gs - STT- LeapSecBDS;  % 신호전달 시간 보정... bds윤초 반영
                end
                SatPos = GetSatPosNC_GC(eph, icol, tc); % 위성위치 산출
                SatPos = RotSatPos(SatPos, STT); % 지구자전효과 고려
                vec_rho = SatPos - vec_site';
                rho = norm(vec_rho);
                [az,el] = xyz2azel(vec_rho, AppLat, AppLon);
                if prn > 100 && prn < 200
                    %                     dIono = 0;
                    %                     dTrop = 0;
                    dIono = ionoKlob(al, be, gs, az, el, x(1:3)); % 이온층 보정(Klobuchar 모델)
                    dTrop = ZHD2SHD(gw, gs, TruePos, el, ZHD); % 대류권 보정
                elseif prn > 200 && prn < 300
                    %                     dIono = 0;
                    %                     dTrop = 0;
                    dIono = ionoKlob(al, be, gs, az, el, x(1:3)); % 이온층 보정(Klobuchar 모델)
                    dTrop = ZHD2SHD(gw, gs, TruePos, el, ZHD); % 대류권 보정
                end
                %             dRel=0;
                dRel = GetRelBRDC(eph, icol, tc); % 상대성효과
                dtSat = a + b*(tc - toe) + c*(tc - toe)^2 - Tgd + dRel; % 위성시계오차 계산... 그룹딜레이, 상대성효과 보정
                if el >=eleCut % 임계고도각
%                     W = 1;
                    W = MakeW_elsnr(el,snr);
                    if prn > 100 && prn < 200
                        %                         PRC = PickPRC(RTCM,prn,gs);
                        PRC = 0;
                        com = rho + x(4) - CCC*dtSat + dTrop + dIono - PRC; % gps 계산값
                        H = [-vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 1 0];
                    elseif prn > 200
                        %                         PRC = PickPRC(RTCM,prn,gs);
                        PRC = 0;
                        com = rho + x(5) - CCC*dtSat + dTrop + dIono - PRC; % bds 계산값
                        H = [-vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 0 1];
                    end
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
                estm(nEst,2:6) = x(1:5);
                estm(nEst,7) = NoSatsUsed;
                Scount(nEst,1) = NoSatsUsed;
%                 fprintf('%8d : %3d : %8.2f : %8.2f\n', j, iter, x(1)' - TruePos(1), x(2)' - TruePos(2));
                break;
            end
        end
    end
    estm = estm(1:nEst,:);
    Scount=Scount(1:nEst, :);
    
    
    fixed_tmp = ([FileQM,'=estm;']);
    eval(fixed_tmp);
    toc;
end

% figure(101);
% [dXYZ, dNEV]=PosErrors(estm(:,1), TruePos,estm(:,2:5),Scount);
% break
%
% %% 입력받은 자료의 크기 및 참값에 해당하는 위경도 결정
% NoPos = length(estm(:,1));
% gd = xyz2gd(TruePos); TrueLat = gd(1); TrueLon = gd(2);
% %% 추정된 XYZ와 참값 XYZ의 차이값 결정
% dXYZ = zeros(NoPos,3);
% for k = 1:NoPos
%     gd(k,:) = xyz2gd(estm(k,2:4));
% end
% estm_L=estm;
% estm_R=estm;
% % toc
% % break
%
%
% open ariport_final2.fig
% hold on
% plot(gd(:,2),gd(:,1),'o:')
% axis([126.45219 126.45255 37.44295 37.44325])
% xlabel('longitude')
% ylabel('latitude')
% title('161223 인천공항 Pole 1')
% break

