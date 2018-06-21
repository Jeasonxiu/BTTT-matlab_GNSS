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
% obsfileBs = '16.obs';                                 % Reference Obsevation
% navfileBs = '16.nav';                                 % Reference Navigation
% obsfileRv = '36.obs';                                 % estm Obsevation
% navfileRv = '36.nav';                                 % estm Navigation
obsfileBs = 'jf092190.obs';                                 % Reference Obsevation
navfileBs = 'jf092190.nav';                                 % Reference Navigation
obsfileRv = 'jr092190.obs';                                 % estm Obsevation
navfileRv = 'jr092190.nav';                                 % estm Navigation
% obsfileBs = 'jfront10.obs';                                 % Reference Obsevation
% navfileBs = 'jfront10.nav';                                 % Reference Navigation
% obsfileRv = 'jrear10.obs';                                 % estm Obsevation
% navfileRv = 'jrear10.nav';                                 % estm Navigation
% obsfileBs = 'DDBs115_33.16o';                                 % Reference Obsevation
% navfileBs = 'DDBs115_33.16n';                                 % Reference Navigation
% obsfileRv = 'DDRv115_13.16o';                                 % estm Obsevation
% navfileRv = 'DDRv115_13.16n';                                 % estm Navigation
% obsfileBs = 'SBBs055_12.16o';                                 % Reference Obsevation
% navfileBs = 'SBBs055_12.16n';                                 % Reference Navigation
% obsfileRv = 'SBRv055_32.16o';                                 % estm Obsevation
% navfileRv = 'SBRv055_32.16n';                                 % estm Navigation
% obsfileBs = 'uf062190g.obs';                                 % Reference Obsevation
% navfileBs = 'uf062190g.nav';                                 % Reference Navigation
% obsfileRv = 'ur062190g.obs';                                 % estm Obsevation
% navfileRv = 'ur062190g.nav';                                 % estm Navigation
%% True Distance
Truedis = 1.41;         % 신항대로 실험
% Truedis = 4.8;
%% 제외 PRN 설정
OPRN = 0;      % 제외 PRN 없는 경우 0
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
QM1 = SelectQM(arrQM1, ObsType); QM1 = QM1(find(QM1(:,2) ~= OPRN),:);
QM11 = SelectQM(arrQM1, ObsType2); QM11 = QM11(find(QM11(:,2) ~= OPRN),:);
QM111 = SelectQM(arrQM1, 131); QM111 = QM111(find(QM111(:,2) ~= OPRN),:);
[arrQM2, FinalPRNs2, FinalTTs2] = ReadQM(renameRv);
QM2 = SelectQM(arrQM2, ObsType); QM2 = QM2(find(QM2(:,2) ~= OPRN),:);
QM22 = SelectQM(arrQM2, ObsType2); QM22 = QM22(find(QM22(:,2) ~= OPRN),:);
QM222 = SelectQM(arrQM2, 131); QM222 = QM222(find(QM222(:,2) ~= OPRN),:);

%% load QM
% load('DDevent.mat');

%% 항법메시지를 읽어들여서 행렬로 저장하고, Klobuchar 모델 추출
eph = ReadEPH(navfileBs);
[al, be] = GetALBE(navfileBs);

%% PRN 제외할 구역 설정
OSTART = 0;
OSTOP = 560445;

%% 라이넥스 파일에서 대략적인 관측소 좌표를 뽑아냄
AppPos = GetAppPos(obsfileRv);
if AppPos(1) == 0
    AppPos = TruePosRv;
end
gd = xyz2gd(AppPos); AppLat = gd(1); AppLon = gd(2);

%% 라이넥스 파일에서 Base Station의 좌표를 계산함
Bs = PP(obsfileBs,navfileBs);              % without Correction
% Bs = PPwC(obsfileBs, navfileBs);            % with Correction
% Bs = PP_prn(obsfileBs, navfileBs, OSTART, OSTOP, OPRN);            % without Correction, 제외 PRN 및 시간
% Bs = PPwC_prn(obsfileBs, navfileBs, OSTART, OSTOP, OPRN);            % with Correction, 제외 PRN 및 시간
Bsxyz = [];
% [Bsgd, Bsgs, Bsutc, Bsla, Bslo, Bsh, Bsxyz] = GGA2gd(strcat(obsfileBs(1:(length(obsfileBs)-4)),'.ubx'));
% Bs = [Bsgs Bsxyz];

%% 선택 시간 sorting
[year, month, days, hour, minute, second]= obs2date(obsfileBs);
[gws, gday] = ydoy2gwgd(YY, DOY);
[gws, gsec] = date2gwgs(year, month, days, hour, minute, second);
% VST = [17, 21, 55];
% VET = [17, 26, 00];
% [gws, start_time] = date2gwgs(year, month, days, VST(1)-9, VST(2), VST(3)); start_time = round(start_time) - 17;
% [gws, end_time] = date2gwgs(year, month, days, VET(1)-9, VET(2), VET(3)); end_time = round(end_time) - 17;
% Bs = Bs(find(Bs(:,1) == start_time):find(Bs(:,1) == end_time),:);

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
    
    gs = FinalTTs(j);
    %% 해당 시각 Bs의 위치 값 찾기
    Base(j,:) = Bs(find(Bs(:,1) == gs),:);
    %     if j == 56
    %         Base(j,2:4) = [-3054757.26302859 	4036638.61111360 	3867124.69787402 ];
    %     end
    %     if j == 57
    %         Base(j,2:4) = [-3054759.50858828 	4036622.00804043 	3867140.64891300 ];
    %     end
    TruePosBs_PP = Base(j,2:4);
    %% base의 gd
    base_gd(j,1) = gs;
    base_gd(j,2:4) = xyz2gd(TruePosBs_PP(1:3)); % Base의 xyz값 gd 로 변환
    AppLatBs = base_gd(j,2); AppLonBs = base_gd(j,3);
    %% 해당 시각 gs의 관측치 선별 및 공통관측 위성 찾기
    indexQM1 = find(QM1(:,1) == gs);
    QM1eBs = QM1(indexQM1,:);       % Base Pseudo-Range(gs)
    QM11eBs = QM11(indexQM1,:);     % Base SNR(gs)
    %% rtklib로 ubx를 변환 시 navfile에 위성정보가 부족할때를 대비하기 위해
    existprnBs = intersect(unique(eph(:,18)), QM1eBs(:,2));
    arrSVBs = zeros(length(existprnBs),1);
    for kk = 1:length(existprnBs)
        arrSVBs(kk) = find(QM1eBs(:,2) == existprnBs(kk,1));
    end
    QM1eBs = QM1eBs(sort(arrSVBs),:);
    QM11eBs = QM11eBs(sort(arrSVBs),:);
    
    indexQM2 = find(QM2(:,1) == gs);
    QM2eRv = QM2(indexQM2,:);           % Rover Pseudo-Range(gs)
    QM22eRv = QM22(indexQM2,:);         % Rover SNR(gs)
    QM222eRv = QM222(indexQM2,:);       % Rover Doppler(gs)
    %% rtklib로 ubx를 변환 시 navfile에 위성정보가 부족할때를 대비하기 위해
    existprnRv = intersect(unique(eph(:,18)), QM2eRv(:,2));
    arrSVRv = zeros(length(existprnRv),1);
    for kkk = 1:length(existprnRv)
        arrSVRv(kkk) = find(QM2eRv(:,2) == existprnRv(kkk,1));
    end
    QM2eRv = QM2eRv(sort(arrSVRv),:);
    QM22eRv = QM22eRv(sort(arrSVRv),:);
    QMdop = QM222eRv(sort(arrSVRv),:);
    
    Sats = intersect(QM1eBs(:, 2), QM2eRv(:, 2));
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
    
    
    %% Base 수신기 상태가 않좋을 경우, t-2,t-1 epoch에서의 X,Y,Z 변화률을 이용하여 t X,Y,Z 예상
    if j > 3 && S1sum(j,1)/length(QM11eBs(:,4)) <= 27 && abs(S1sum(j,1)-S1sum(j,2)) >= 100
        %         Base(j,2:4) = Base(j-1,2:4) + (Base(j-1,2:4)-Base(j-3,2:4))/3;
        Base(j,2:4) = Base(j-1,2:4) - Bs(j,6:8);
        TruePosBs_PP = Base(j,2:4);
    elseif j > 3 && S1sum(j,1)/length(QM11eBs(:,4)) <= 32
        %         Base(j,2:4) = Base(j-1,2:4) + (Base(j-1,2:4)-Base(j-3,2:4))/3;
        Base(j,2:4) = Base(j-1,2:4) - Bs(j,6:8);
        TruePosBs_PP = Base(j,2:4);
    elseif j > 3 && abs(S1sum(j,1)-S1sum(j,2)) >= 100
        %         Base(j,2:4) = Base(j-1,2:4) + (Base(j-1,2:4)-Base(j-3,2:4))/3;
        Base(j,2:4) = Base(j-1,2:4) - Bs(j,6:8);
        TruePosBs_PP = Base(j,2:4);
    end
    
    
    Sats = intersect(QM1eBs(:, 2), QM2eRv(:, 2));
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
            
%             W =1;
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
            estm(nEst,5) = NoSats;
            estm(nEst,6) = NoSatsUsed;  % : 기준위성을 포함해야 함
            estm(nEst,7) = 0;  % : snr issue 가 없을때
            
            
            
            break;
        end
        
    end
    %     QMdop = QMdop(indxUsedSat,:);
    %     estm(nEst,8:10) = GetVelDop(gs, QMdop(1:NoSatsUsed,2), QMdop(1:NoSatsUsed,4), eph, x(1:3)'); % doppler vel
    estm(nEst,8:10) = GetVelDop(gs, QMdop(1:NoSats,2), QMdop(1:NoSats,4), eph, x(1:3)'); % doppler vel
    S1sum(nEst,:) = S1sum(j,:);
    base(nEst,:) = Base(j,:);
    %% snr isuue 발생시 rover 위치 좌표 예측
    if j > 3 && S1sum(nEst,1)/length(QM11eBs(:,4)) <= 32 && abs(S1sum(nEst,1)-S1sum(nEst,2)) >= 100
        %         estm(nEst,2:4) = estm(nEst-1,2:4) + (estm(nEst-1,2:4)-estm(nEst-3,2:4))/3;      % : rover의 XYZ 변화량
        %         estm(nEst,2:4) = estm(nEst-1,2:4) + (Base(nEst-1,2:4)-Base(nEst-3,2:4))/3;      % : base의 XYZ 변화량
        if estm(nEst,1) - estm(nEst-1,1) <= 2
            %             estm(nEst,2:4) = estm(nEst-1,2:4) - estm(nEst,8:10);       % : rover의 도플러 속도 XYZ 변화량
            estm(nEst,2:4) = estm(nEst-1,2:4) - base(nEst,6:8);            % : base의 도플러 속도 XYZ 변화량
            estm(nEst,7) = 1;end
    elseif j > 3 && S1sum(nEst,1,1)/length(QM11eBs(:,4)) <= 32
        %             estm(nEst,2:4) = estm(nEst-1,2:4) + (estm(nEst-1,2:4)-estm(nEst-3,2:4))/3;      % : rover의 XYZ 변화량
        %             estm(nEst,2:4) = estm(nEst-1,2:4) + (Base(nEst-1,2:4)-Base(nEst-3,2:4))/3;      % : base의 XYZ 변화량
        if estm(nEst,1) - estm(nEst-1,1) <= 2
            %             estm(nEst,2:4) = estm(nEst-1,2:4) - estm(nEst,8:10);       % : rover의 도플러 속도 XYZ 변화량
            estm(nEst,2:4) = estm(nEst-1,2:4) - base(nEst,6:8);            % : base의 도플러 속도 XYZ 변화량
            estm(nEst,7) = 2;end
    elseif j > 3 && abs(S1sum(nEst,1)-S1sum(nEst,2)) >= 100
        %             estm(nEst,2:4) = estm(nEst-1,2:4) + (estm(nEst-1,2:4)-estm(nEst-3,2:4))/3;      % : rover의 XYZ 변화량
        %             estm(nEst,2:4) = estm(nEst-1,2:4) + (Base(nEst-1,2:4)-Base(nEst-3,2:4))/3;      % : base의 XYZ 변화량
        if estm(nEst,1) - estm(nEst-1,1) <= 2
            %             estm(nEst,2:4) = estm(nEst-1,2:4) - estm(nEst,8:10);       % : rover의 도플러 속도 XYZ 변화량
            estm(nEst,2:4) = estm(nEst-1,2:4) - base(nEst,6:8);            % : base의 도플러 속도 XYZ 변화량
            estm(nEst,7) = 3;end
    end
    %% rover's longi, lati
    rover_gd(nEst,1) = gs;
    rover_gd(nEst,2:4) = xyz2gd(estm(nEst,2:4)); % rover의 xyz값 gd 로 변환
    AppLat = rover_gd(nEst,2); AppLon = rover_gd(nEst,3);
    
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
[DDdXYZ, DDdNEV, DDdis, DDrms, DDstd] = PostErrorsDDmv(estm, Base, Truedis,0, 10, S1sum);         % 임의의 장소에서 이동 측위시
% [QMnewBs, QMnewRv, QMnew] = DDSkyplot(QM1, QM2, eph, Base, estm);               % Skyplot
% DDPlotQM(renameBs, renameRv, 141)

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

%% 각 위성별 SNR Plot
% DDPlotQM(QM11, QM22, 141, 'Base', 'Rover')








