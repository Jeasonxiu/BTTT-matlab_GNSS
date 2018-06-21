%% 코드의사거리 이중차분 알고리즘
% 07/01/2016 : Joonseong
close all; clear all;

%% 불변 변수 설정 : 빛의 속도, 관측치
CCC = 299792458.;   % CCC = Speed of Light [m/s]
ObsType = 320;      % 사용할 관측지 설정 - 120 : C/A = C1
ObsType2 = 341;      % 사용할 관측치 설정 - 141: snr = S1

%% 임계고도각 설정
eleCut = 15;

%% YY, DOY 설정
YY = 16;
DOY = 50;

%% True Distance
Truedis = 1.41;

%% read SP3
file_sp3 = 'iac19073.sp3';
% file_sp3 = 'igl18845.sp3';    % SBS 강화 vs JPRT
[arrSP3] = ReadSP3_GLO(file_sp3);

%% QMfile load
renameBs = 'QMjfr2';
renameRv = 'QMjrr2';
% renameBs = 'QSBSA_16050a';      % SBS 강화 vs JPRT  
% renameRv = 'QJPRT_16050a';      % SBS 강화 vs JPRT
%% 기준위성 및 제외 PRN 설정
RS = 23;
OPRN = 0;      % 제외 PRN 없는 경우 0

%% QM 파일 읽어들여서 행렬로 저장하고, 사용할 관측치 추출
[arrQM1, FinalPRNs1, FinalTTs1] = ReadQM(renameBs);
QM1 = SelectQM(arrQM1, ObsType); QM1 = QM1(find(QM1(:,2) ~= OPRN),:);
QM11 = SelectQM(arrQM1, ObsType2); QM11 = QM11(find(QM11(:,2) ~= OPRN),:);
QM111 = SelectQM(arrQM1, 331); QM111 = QM111(find(QM111(:,2) ~= OPRN),:);
[arrQM2, FinalPRNs2, FinalTTs2] = ReadQM(renameRv);
QM2 = SelectQM(arrQM2, ObsType); QM2 = QM2(find(QM2(:,2) ~= OPRN),:);
QM22 = SelectQM(arrQM2, ObsType2); QM22 = QM22(find(QM22(:,2) ~= OPRN),:);
QM222 = SelectQM(arrQM2, 331); QM222 = QM222(find(QM222(:,2) ~= OPRN),:);

%% Base 좌표 계산
obsfileBs = 'jf022090r.obs';                                 % Reference Obsevation
navfileBs = 'jf022090r.nav';                                 % Reference Navigation
Bs = PP(obsfileBs,navfileBs);              % without Correction
AppPos = GetAppPos(obsfileBs);
% AppPos = [-3026795.499 4067267.161 3857084.459];
Bs = PP_glo(file_sp3, AppPos, renameBs);            % with Correction

% load('gloBs.mat');
%% PRN 제외할 구역 설정
OSTART = 0;
OSTOP = 560445;

%% 선택 시간 sorting
% [year, month, days, hour, minute, second]= obs2date(obsfileBs);
[gws, gday] = ydoy2gwgd(YY, DOY);
% [gws, gsec] = date2gwgs(year, month, days, hour, minute, second);
% VST = [17, 21, 55];
% VET = [17, 26, 00];
% [gws, start_time] = date2gwgs(year, month, days, VST(1)-9, VST(2), VST(3)); start_time = round(start_time) - 17;
% [gws, end_time] = date2gwgs(year, month, days, VET(1)-9, VET(2), VET(3)); end_time = round(end_time) - 17;
% Bs = Bs(find(Bs(:,1) == start_time):find(Bs(:,1) == end_time),:);

%% 두 QM 파일에서 공통시각(epoch) 추출
FinalTTs = intersect(Bs(:, 1), QM2(:, 1));

%% SBS 강화 vs JPRT 실행시
% FinalTTs = intersect(QM1(:, 1), QM2(:, 1));   % SBS 강화 vs JPRT
% Bs = zeros(3600,9); Bs(:,1) = FinalTTs;       % SBS 강화 vs JPRT

%% 라이넥스 파일에서 대략적인 관측소 좌표를 뽑아냄
AppPos = GetAppPos(obsfileBs);
if AppPos(1) == 0
    AppPos = TruePosRv;
end
gd = xyz2gd(AppPos); AppLat = gd(1); AppLon = gd(2);

%% 추정에 필요한 초기치 설정
MaxIter = 4;
EpsStop = 1e-6;
x = AppPos';


%% 추정과정 시작
NoEpochs = length(FinalTTs);
estm = zeros(NoEpochs, 6);  % : c1(gs), c2/3/4(X/Y/Z of Rv), c5(#Sats_obs), c6(#Sats_used)
nEst = 0;
cnt = 0;
No_Sat = 0;
for j = 1:NoEpochs
%     for j = 1:600
%     for j = 130:260
    
    gs = FinalTTs(j);
    %% 해당 시각 Bs의 위치 값 찾기
    Base(j,:) = Bs(find(Bs(:,1) == gs),:);
    Base(j,2:4) = [-3003051.372 4059906.090 3883100.144];       %% SBS 강화 좌표
    TruePosBs_PP = Base(j,2:4);

    %% base의 gd
    base_gd(j,1) = gs;
    base_gd(j,2:4) = xyz2gd(TruePosBs_PP(1:3)); % Base의 xyz값 gd 로 변환
    AppLatBs = base_gd(j,2); AppLonBs = base_gd(j,3);
    %% 해당 시각 gs의 관측치 선별 및 공통관측 위성 찾기
    indexQM1 = find(QM1(:,1) == gs);
    QM1eBs = QM1(indexQM1,:);       % Base Pseudo-Range(gs)
    QM11eBs = QM11(indexQM1,:);     % Base SNR(gs)
       
    indexQM2 = find(QM2(:,1) == gs);
    QM2eRv = QM2(indexQM2,:);           % Rover Pseudo-Range(gs)
    QM22eRv = QM22(indexQM2,:);         % Rover SNR(gs)
    QMdop = QM222(indexQM2,:);       % Rover Doppler(gs)
    
    Sats = intersect(QM1eBs(:, 2), QM2eRv(:, 2));
%     for k = 1: length(Sats)
%         arrglosat(k,1:4) = IntpSP3e1(arrSP3, Sats(k), gs);
%         arrglosat(k,5) = Sats(k,1);
%     end
%% 기준위성 RS와 다른위성 OS 선별/ SatsEl - c1(gs), c2(prn), c3(el)
    [SatsEl, indxRS] = PickRSelGloSP3(gs, Sats, arrSP3, TruePosBs_PP);  % : RS 기준위성
    RS = Sats(indxRS); RefSV(j,1) = RS;

    NoSats = length(Sats); No_Sat = No_Sat + NoSats;
    if OSTART == 0
    elseif (gs >= OSTART) && (gs <= OSTOP)
        QM1eBs = QM1eBs(find(QM11eBs(:,2) ~= OPRN),:);
        QM2eRv = QM2eRv(find(QM22eRv(:,2) ~= OPRN),:);
        QM11eBs = QM11eBs(find(QM11eBs(:,2) ~= OPRN),:);
        QM22eRv = QM22eRv(find(QM22eRv(:,2) ~= OPRN),:);
    end
    
    %% Base 수신 위성의 전체 SNR 합을 구하고 평균값이용
    S1sum(j,1) = sum(QM11eBs(:,4));
    S1sum(j,2) = sum(QM22eRv(:,4));
    
    if j > 3 && S1sum(j,1)/length(QM11eBs(:,4)) <= 32 && abs(S1sum(j,1)-S1sum(j,2)) >= 100
        S1sum(j,3) = 1;
    elseif j > 3 && S1sum(j,1)/length(QM11eBs(:,4)) <= 32
        S1sum(j,3) = 2;
    elseif j > 3 && abs(S1sum(j,1)-S1sum(j,2)) >= 100
        S1sum(j,3) = 3;
    else
        S1sum(j,3) = 0;
    end
    S1sum(j,4) = length(QM11eBs(:,4));
    S1sum(j,5) = length(QM22eRv(:,4));
    
    %% 기준위성 좌표 먼저 계산 - Bs 기준
    STT = GetSTTsp3(gs, RS, arrSP3, TruePosBs_PP');
    tc = gs - STT;
    vec_RS = IntpSP3e1(arrSP3, RS, tc); vec_RS = vec_RS(2:4);
    S1BsRS = QM11eBs(find(QM11eBs(:,2) == RS), 4);      % BASE RS SNR matrix
    S1RvRS = QM22eRv(find(QM22eRv(:,2) == RS), 4);      % ROVER RS SNR matrix
    
    for Iter = 1:MaxIter
        
        HTH = zeros(3,3);
        HTy = zeros(3,1);
        NoSatsUsed = NoSats;
        usedSatCount = 0;
        %% 각 위성에 대한 관측치, 계산치, H행렬 계산
        for kS = 1:NoSats
            if kS == indxRS || SatsEl(kS, 3) < eleCut
                if SatsEl(kS, 3) < eleCut
                    NoSatsUsed = NoSatsUsed - 1;
                    disp([RS SatsEl(kS, 2)])
                end
                continue
            end
            cnt = cnt + 1;
            %% DD 관측치 생성 파트 --- 추후에 for 루프 밖으로 빼야 함 11/8/14
            OS = Sats(kS);
            if Iter == 1
                if OS > 0
                    usedSatCount = usedSatCount + 1;
                    indxUsedSat(usedSatCount,1) = find(QMdop(:,2) == OS);
                    if NoSatsUsed - usedSatCount == 1
                        indxUsedSat(usedSatCount+1,1) = find(QMdop(:,2) == RS);
                    end
                end
            end
            OtherSats(kS,1) = OS;
            
            S1BsOS = QM11eBs(find(QM11eBs(:,2) == OS), 4);      % BASE OS SNR matrix
            S1RvOS = QM22eRv(find(QM22eRv(:,2) == OS), 4);      % ROVER OS SNR matrix
            
            obs_BsRS = QM1eBs(find(QM1eBs(:, 2) == RS), 4);
            obs_RvRS = QM2eRv(find(QM2eRv(:, 2) == RS), 4);
            obs_BsOS = QM1eBs(find(QM1eBs(:, 2) == OS), 4);
            obs_RvOS = QM2eRv(find(QM2eRv(:, 2) == OS), 4);
            obs = (obs_BsRS - obs_RvRS) - (obs_BsOS - obs_RvOS);
            %% DD 계산치 생성 파트 - 기타위성 좌표 계산(기준위성 좌표는 이미 계산 완료)
            STT = GetSTTsp3(gs, OS, arrSP3, x(1:3));
            tc = gs - STT;
            vec_OS = IntpSP3e1(arrSP3, OS, tc); vec_OS = vec_OS(2:4);
            
            %% DD 계산치 생성 파트 - 거리계산치를 각각 계산한 다음 DD 계산치 계산
            vec_BsRS = vec_RS - TruePosBs_PP;  com_BsRS = norm(vec_BsRS);
            vec_RvRS = vec_RS - x(1:3)';    com_RvRS = norm(vec_RvRS);
            vec_BsOS = vec_OS - TruePosBs_PP;  com_BsOS = norm(vec_BsOS);
            vec_RvOS = vec_OS - x(1:3)';    com_RvOS = norm(vec_RvOS);
            com = (com_BsRS - com_RvRS) - (com_BsOS - com_RvOS);
            y = obs -com;
            %% 각 위성 az, el 저장
            
            [azBsRS,elBsRS] = xyz2azel(vec_BsRS, AppLatBs, AppLonBs);
            azelBsRS(j,:) = [azBsRS, elBsRS];
            [azBsOS,elBsOS] = xyz2azel(vec_BsOS, AppLatBs, AppLonBs);
            azelBsOS(cnt,:) = [azBsOS, elBsOS];
            %% Weighting
            S1RS(cnt,:) = [S1BsRS, S1RvRS,gs,OS];
            S1OS(cnt,:) = [S1BsOS, S1RvOS,gs,OS];
            DDel(cnt,:) = [elBsRS, elBsOS,gs,OS];

            W =1;
            W = DDMakeW_elsnr(S1RS(cnt,:),S1OS(cnt,:),DDel(cnt,:));
            weight(cnt,:) = W;
            
            %% H 행렬 계산 파트
            H(1,1) = vec_RvRS(1)/com_RvRS - vec_RvOS(1)/com_RvOS;
            H(1,2) = vec_RvRS(2)/com_RvRS - vec_RvOS(2)/com_RvOS;
            H(1,3) = vec_RvRS(3)/com_RvRS - vec_RvOS(3)/com_RvOS;
            
            HTH = HTH + H'*W*H;
            HTy = HTy + H'*W*y;
        end
        OTHER{j,1} = OtherSats;
        OTHER{j,2} = indxUsedSat;
        xhat = inv(HTH) * HTy;
        x = x + xhat;
        
        if norm(xhat) < EpsStop;
            nEst = nEst + 1;
            estm(nEst,1) =gs;
            estm(nEst,2:4) =x;
%             estm(nEst,2:4) =x + [23744.127 	-7361.071 	26015.685]' ;       % SBS 강화 vs JPRT
            estm(nEst,5) = NoSats;
            estm(nEst,6) = NoSatsUsed;  % : 기준위성을 포함해야 함
            estm(nEst,7) = 0;  % : snr issue 가 없을때
            
            
            
            break;
        end
        
    end
    
    S1sum(nEst,:) = S1sum(j,:);
    base(nEst,:) = Base(j,:);
    %% snr isuue 발생시 rover 위치 좌표 예측
%     if j > 3 && S1sum(nEst,1)/length(QM11eBs(:,4)) <= 32 && abs(S1sum(nEst,1)-S1sum(nEst,2)) >= 100
%         %         estm(nEst,2:4) = estm(nEst-1,2:4) + (estm(nEst-1,2:4)-estm(nEst-3,2:4))/3;      % : rover의 XYZ 변화량
%         %         estm(nEst,2:4) = estm(nEst-1,2:4) + (Base(nEst-1,2:4)-Base(nEst-3,2:4))/3;      % : base의 XYZ 변화량
%         if estm(nEst,1) - estm(nEst-1,1) <= 2
%             %             estm(nEst,2:4) = estm(nEst-1,2:4) - estm(nEst,8:10);       % : rover의 도플러 속도 XYZ 변화량
%             estm(nEst,2:4) = estm(nEst-1,2:4) - base(nEst,6:8);            % : base의 도플러 속도 XYZ 변화량
%             estm(nEst,7) = 1;end
%     elseif j > 3 && S1sum(nEst,1,1)/length(QM11eBs(:,4)) <= 32
%         %             estm(nEst,2:4) = estm(nEst-1,2:4) + (estm(nEst-1,2:4)-estm(nEst-3,2:4))/3;      % : rover의 XYZ 변화량
%         %             estm(nEst,2:4) = estm(nEst-1,2:4) + (Base(nEst-1,2:4)-Base(nEst-3,2:4))/3;      % : base의 XYZ 변화량
%         if estm(nEst,1) - estm(nEst-1,1) <= 2
%             %             estm(nEst,2:4) = estm(nEst-1,2:4) - estm(nEst,8:10);       % : rover의 도플러 속도 XYZ 변화량
%             estm(nEst,2:4) = estm(nEst-1,2:4) - base(nEst,6:8);            % : base의 도플러 속도 XYZ 변화량
%             estm(nEst,7) = 2;end
%     elseif j > 3 && abs(S1sum(nEst,1)-S1sum(nEst,2)) >= 100
%         %             estm(nEst,2:4) = estm(nEst-1,2:4) + (estm(nEst-1,2:4)-estm(nEst-3,2:4))/3;      % : rover의 XYZ 변화량
%         %             estm(nEst,2:4) = estm(nEst-1,2:4) + (Base(nEst-1,2:4)-Base(nEst-3,2:4))/3;      % : base의 XYZ 변화량
%         if estm(nEst,1) - estm(nEst-1,1) <= 2
%             %             estm(nEst,2:4) = estm(nEst-1,2:4) - estm(nEst,8:10);       % : rover의 도플러 속도 XYZ 변화량
%             estm(nEst,2:4) = estm(nEst-1,2:4) - base(nEst,6:8);            % : base의 도플러 속도 XYZ 변화량
%             estm(nEst,7) = 3;end
%     end
    %% rover's longi, lati
%     rover_gd(nEst,1) = gs;
%     rover_gd(nEst,2:4) = xyz2gd(estm(nEst,2:4)); % rover의 xyz값 gd 로 변환
%     AppLat = rover_gd(nEst,2); AppLon = rover_gd(nEst,3);
    
end

%% 측위오차 분석 & 그래프 작성
estm = estm(1:nEst, :);
% estm(56,2:4) = [-3054756.96353402 	4036641.42409425 	3867120.42573697 ];
% estm(57,2:4) = [-3054758.84268655 	4036624.69832510 	3867136.75460206  ];
% estm(52:58,2:4) = [-3054752.00611902 	4036707.92645826 	3867056.26333367;...
%     -3054753.38719373 	4036690.37587079 	3867071.21295368;...
%     -3054753.43609467 	4036674.29539320 	3867088.48809783;...
%     -3054755.08438149 	4036658.14986339 	3867104.09687187 ;...
%     -3054756.96353402 	4036641.42409425 	3867120.42573697 ;...
%     -3054758.84268655 	4036624.69832510 	3867136.75460206 ;...
%     -3054760.72183908 	4036607.97255596 	3867153.08346716];

% [dXYZ, dNEV] = PosTErrorsJOON(estm(:, 1), TruePosRv, estm(:, 2:5));         % A,B Point 정지 측위시
% [DDdXYZ, DDdNEV, DDdis, DDrms, DDstd] = PostErrorsDD(estm, Base, TruePosBs, TruePosRv); % A,B Point 정지 측위시
% [DDdXYZ, DDdNEV, DDdis, DDrms, DDstd] = PostErrorsDDmv(estm, Base, Truedis,0, 5, S1sum);         % 임의의 장소에서 이동 측위시
% [QMnewBs, QMnewRv, QMnew] = DDSkyplot(QM1, QM2, eph, Base, estm);               % Skyplot
% DDPlotQM(renameBs, renameRv, 141)

%% SBS 강화 vs JPRT 실행시
[dXYZ, dNEV] = PosTErrors2(estm(:,1),  [-3026795.499 4067267.161 3857084.459], estm(:,2:5));

figure(99)
subplot(3,2,5)
hold on; grid on;
xlim([0 length(estm)])
plot(S1sum(:,1),'r.:');
plot(S1sum(:,2),'b.:');
legend('Base S1 sum','Rover S1 sum')



%% 특이 지점 시간 탐색
for TT = 1:length(estm(:,1))
    event = DDdis(TT,2);
    if event - Truedis > 1.5
        eT = estm(TT,1)-17;
        [yyy, mmo, ddd, hhh, mmm, sss] = gwgs2date(gws, eT);
        event_time(TT,:) = [hhh+9,mmm,sss,event-Truedis,estm(TT,1),1];
    else
        eT = estm(TT,1)-17;
        [yyy, mmo, ddd, hhh, mmm, sss] = gwgs2date(gws, eT);
        event_time(TT,:) = [hhh+9,mmm,sss,event-Truedis,estm(TT,1),0];
    end
end

% PPlotQM(QM11, QM22, 141, 'Base', 'Rover')
%% 동일 epoch base, rover 구글 Plot
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
% %     plot(Bs_gd(2), Bs_gd(1),'r.','MarkerSize',20)
% %     plot(Rv_gd(2), Rv_gd(1),'b.','MarkerSize',20)
% %     legend('epoch','Left','Right')
% end
%     figure(200)
%     grid on; hold on;
%     plot(base_gd(:,3), base_gd(:,2),'r.','MarkerSize',10)
%     plot(rover_gd(:,3), rover_gd(:,2),'b.','MarkerSize',10)
% DDPlotQMcheck(renameBs, renameRv, 141,min(Bs(:,1)),max(Bs(:,1)),event_time)

% for i = 1:length(DDdis)
%     DDdis2D(i,:) = DDdis(i,1);
%     DDdis3D(i,:) = DDdis(i,2);
%     TrueDis(i,:) = DDdis(i,3);
%     figure(999)
%     hold on; grid on;
%     xlim([0 length(estm)])
%     plot(DDdis2D,'r-');
% %     plot(DDdis3D,'b-');
%     plot(TrueDis,'k-');
%     drawnow
% end
%
% figure(999)
% hold on; grid on;
% xlabel({['Double Differencing dNE = ', num2str(DDrms(1)), '   std =', num2str(DDstd(1))]});
% % xlabel({['Double Differencing dNE = ', num2str(DDrms(1)), '   std =', num2str(DDstd(1))],...
% %     ['Double Differencing  3D = ', num2str(DDrms(2)), '   std =', num2str(DDstd(2))]});
% ylabel('Distance(meter)');
% legend('2D(dNE) Distance')
% % legend('2D(dNE) Distance','3D(d3D) Distance')