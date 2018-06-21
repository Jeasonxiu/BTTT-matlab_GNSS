%% 코드의사거리 이중차분 알고리즘
% 07/01/2016 : Joonseong
close all; clear all;

%% 불변 변수 설정 : 빛의 속도, 관측치
CCC = 299792458.;   % CCC = Speed of Light [m/s]
ObsType = 120;      % 사용할 관측지 설정 - 120 : C/A = C1
ObsType2 = 141;      % 사용할 관측치 설정 - 141: snr = S1

%% 임계고도각 설정
eleCut = 15;

%% QM 파일 핸들링
% obsfileBs = 'DDBs115_31.16o';                                 % Reference Obsevation
% navfileBs = 'DDBs115_31.16n';                                 % Reference Navigation
% obsfileRv = 'DDRv115_11.16o';                                 % estm Obsevation
% navfileRv = 'DDRv115_11.16n';                                 % estm Navigation
obsfileBs = 'SBBs055_14.16o';                                 % Reference Obsevation
navfileBs = 'SBBs055_14.16n';                                 % Reference Navigation
obsfileRv = 'SBRv055_34.16o';                                 % estm Obsevation
navfileRv = 'SBRv055_34.16n';                                 % estm Navigation

%% Rinex to QM
WriteObs(obsfileBs);                                        % Base Station QMfile                 
renameBs = renameQMfile(obsfileBs);                    % Base Station's QMfile Renaming
WriteObs(obsfileRv);                                         % estm QMfile
renameRv = renameQMfile(obsfileRv);                    % estm's QMfile Renaming

%% 항법 RINEX 파일 설정

%% 기타 설정 : 사이트 좌표 참값
% TruePosBs = [];                                       % : Not exist True position
% TruePosRv = [];                                       % : Not exist True position
TruePosBs = [-3041235.578 4053941.677 3859881.013];   % : JPspace A point
TruePosRv = [-3041241.741 4053944.143 3859873.640];   % : JPspace B point
[YY, DOY] = obs2YYDOY(obsfileBs);

%% QM 파일 읽어들여서 행렬로 저장하고, 사용할 관측치 추출
[arrQM1, FinalPRNs1, FinalTTs1] = ReadQM(renameBs);
QM1 = SelectQM(arrQM1, ObsType);
QM11 = SelectQM(arrQM1, ObsType2);
[arrQM2, FinalPRNs2, FinalTTs2] = ReadQM(renameRv);
QM2 = SelectQM(arrQM2, ObsType);
QM22 = SelectQM(arrQM2, ObsType2);

%% load QM
% load('DDevent.mat');

%% 항법메시지를 읽어들여서 행렬로 저장하고, Klobuchar 모델 추출
eph = ReadEPH(navfileBs);
[al, be] = GetALBE(navfileBs);

%% 라이넥스 파일에서 대략적인 관측소 좌표를 뽑아냄
AppPos = GetAppPos(obsfileRv);
if AppPos(1) == 0
    AppPos = TruePosRv;
end
gd = xyz2gd(AppPos); AppLat = gd(1); AppLon = gd(2);

%% 라이넥스 파일에서 Base Station의 좌표를 계산함
Bs = PP(obsfileBs,navfileBs);              % without Correction
% Bs = PPwC(obsfileBs, navfileBs);            % with Correction
Bsxyz = [];
% [Bsgd, Bsgs, Bsutc, Bsla, Bslo, Bsh, Bsxyz] = GGA2gd(strcat(obsfileBs(1:(length(obsfileBs)-4)),'.ubx'));
% Bs = [Bsgs Bsxyz];

%% 선택 시간 sorting
% [year, month, days]= obs2date(obsfileBs);
% VST = [17, 12, 25];
% VET = [17, 14, 05];
% [gws, start_time] = date2gwgs(year, month, days, VST(1)-9, VST(2), VST(3)); start_time = round(start_time) - 17;
% [gws, end_time] = date2gwgs(year, month, days, VET(1)-9, VET(2), VET(3)); end_time = round(end_time) - 17;
% Bs = Bs(find(Bs(:,1) == start_time):find(Bs(:,1) == end_time),:);
%% True Distance
Truedis = 3.8;
%% 두 QM 파일에서 공통시각(epoch) 추출
if ~isempty(Bsxyz)
    FinalTTs = intersect(Bsgs(:, 1), QM2(:, 1));
else
    FinalTTs = intersect(Bs(:, 1), QM2(:, 1));
end

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
% for j = 28:28
    
    gs = FinalTTs(j);
    %% 해당 시각 Bs의 위치 값 찾기
    Base(j,:) = Bs(find(Bs(:,1) == gs),1:4);
    TruePosBs_PP = Base(j,2:4);
    %% base의 gd
    base_gd(j,:) = xyz2gd(TruePosBs_PP(1:3)); % Base의 xyz값 gd 로 변환
    AppLatBs = base_gd(j,1); AppLonBs = base_gd(j,2);
    %% 해당 시각 gs의 관측치 선별 및 공통관측 위성 찾기
    indexQM1 = find(QM1(:,1) == gs);
    QM1eBs = QM1(indexQM1,:);
    QM11eBs = QM11(indexQM1,:);
    %% rtklib로 ubx를 변환 시 navfile에 위성정보가 부족할때를 대비하기 위해 
    existprnBs = intersect(unique(eph(:,18)), QM1eBs(:,2));
    arrSVBs = zeros(length(existprnBs),1);
    for kk = 1:length(existprnBs)
        arrSVBs(kk) = find(QM1eBs(:,2) == existprnBs(kk,1));
    end
    QM1eBs = QM1eBs(sort(arrSVBs),:);
    QM11eBs = QM11eBs(sort(arrSVBs),:);
    
    indexQM2 = find(QM2(:,1) == gs);
    QM1eRv = QM2(indexQM2,:);
    QM11eRv = QM22(indexQM2,:);
    %% rtklib로 ubx를 변환 시 navfile에 위성정보가 부족할때를 대비하기 위해
    existprnRv = intersect(unique(eph(:,18)), QM1eRv(:,2));
    arrSVRv = zeros(length(existprnRv),1);
    for kkk = 1:length(existprnRv)
        arrSVRv(kkk) = find(QM1eRv(:,2) == existprnRv(kkk,1));
    end
    QM1eRv = QM1eRv(sort(arrSVRv),:);
    QM11eRv = QM11eRv(sort(arrSVRv),:);
    
    
    Sats = intersect(QM1eBs(:, 2), QM1eRv(:, 2));
    NoSats = length(Sats); No_Sat = No_Sat + NoSats;
    
    %% 기준위성 RS와 다른위성 OS 선별/ SatsEl - c1(gs), c2(prn), c3(el)
    [SatsEl, indxRS] = PickRSel(gs, Sats, eph, TruePosBs_PP);  % : RS 기준위성
    RS = Sats(indxRS); RefSV(j,1) = RS;
    
    %% 기준위성 좌표 먼저 계산 - Bs 기준
    icol = PickEPH(eph, RS ,gs);
    STT = GetSTTbrdc(gs, RS, eph, TruePosBs_PP);
    tc = gs -STT;
    vec_RS = GetSatPosNC(eph, icol, tc);
    vec_RS = RotSatPos(vec_RS, STT);
    S1BsRS = QM11eBs(find(QM11eBs(:,2) == RS), 4);
    S1RvRS = QM11eRv(find(QM11eRv(:,2) == RS), 4);
    for Iter = 1:MaxIter
        
        HTH = zeros(3,3);
        HTy = zeros(3,1);
        NoSatsUsed = NoSats;
        
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
            S1BsOS = QM11eBs(find(QM11eBs(:,2) == OS), 4);
            S1RvOS = QM11eRv(find(QM11eRv(:,2) == OS), 4);
            obs_BsRS = QM1eBs(find(QM1eBs(:, 2) == RS), 4);
            obs_RvRS = QM1eRv(find(QM1eRv(:, 2) == RS), 4);
            obs_BsOS = QM1eBs(find(QM1eBs(:, 2) == OS), 4);
            obs_RvOS = QM1eRv(find(QM1eRv(:, 2) == OS), 4);
            obs = (obs_BsRS - obs_RvRS) - (obs_BsOS - obs_RvOS);
            %% DD 계산치 생성 파트 - 기타위성 좌표 계산(기준위성 좌표는 이미 계산 완료)
            icol = PickEPH(eph, OS, gs);
            STT = GetSTTbrdc(gs, OS, eph, x(1:3)'); % :  OS 위성위치를 추정좌표 기준으로 변경 11/9/14
            tc = gs - STT;
            vec_OS = GetSatPosNC(eph, icol, tc);
            vec_OS = RotSatPos(vec_OS, STT);
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
            W = DDMakeW_elsnr(S1RS(cnt,:),S1OS(cnt,:),DDel(cnt,:));
            weight(cnt,:) = W;
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
            estm(nEst,5) = NoSats;
            estm(nEst,6) = NoSatsUsed;  % : 기준위성을 포함해야 함
            break;
        end
    end
    rover_gd(j,:) = xyz2gd(estm(j,2:4)); % rover의 xyz값 gd 로 변환
    AppLat = rover_gd(j,1); AppLon = rover_gd(j,2);

end

%% 측위오차 분석 & 그래프 작성
estm = estm(1:nEst, :);

% [dXYZ, dNEV] = PosTErrorsJOON(estm(:, 1), TruePosRv, estm(:, 2:5));         % A,B Point 정지 측위시
% [DDdXYZ, DDdNEV, DDdis, DDrms, DDstd] = PostErrorsDD(estm, Base, TruePosBs, TruePosRv); % A,B Point 정지 측위시
[DDdXYZ, DDdNEV, DDdis, DDrms, DDstd] = PostErrorsDDmv(estm, Base, Truedis,0,10);         % 임의의 장소에서 이동 측위시  
% [QMnewBs, QMnewRv, QMnew] = DDSkyplot(QM1, QM2, eph, Base, estm);               % Skyplot
% DDPlotQM(renameBs, renameRv, 141)

%% 특이 지점 시간 탐색
for TT = 1:length(estm(:,1))
    event = DDdis(TT,2);
    if event - Truedis > 1.5
        eT = estm(TT,1)-17;
        [yyy, mmo, ddd, hhh, mmm, sss] = gwgs2date(gws, eT);
        event_time(TT,:) = [hhh+9,mmm,sss,event-Truedis,1];
    else
        eT = estm(TT,1)-17;
        [yyy, mmo, ddd, hhh, mmm, sss] = gwgs2date(gws, eT);
        event_time(TT,:) = [hhh+9,mmm,sss,event-Truedis,0];
    end
end

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