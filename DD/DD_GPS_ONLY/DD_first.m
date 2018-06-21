%% 코드의사거리 이중차분 알고리즘
% 07/01/2016 : Joonseong

%% 불변 변수 설정 : 빛의 속도, 관측치
CCC = 299792458.;   % CCC = Speed of Light [m/s]
ObsType = 120;      % 사용할 관측지 설정 - 120 : C/A = C1

%% 임계고도각 설정
eleCut = 15;

%% QM 파일 핸들링
FileQM1 = 'QZBLA_15334';
FileQM2 = 'QZBLB_15334';

%% 항법 RINEX 파일 설정
FileNav = 'brdc3340.15n';

%% 기타 설정 : 사이트 좌표 참값 & GPS Week
TruePosBs = [-3026789.236 4067255.523 3857098.106]; % : 15334 참값
TruePosRv = [-3026789.236 4067255.523 3857098.106]; % : 15334 참값

%% QM 파일 읽어들여서 행렬로 저장하고, 사용할 관측치 추출
[arrQM1, FinalPRNs1, FinalTTs1] = ReadQM(FileQM1);
QM1 = SelectQM(arrQM1, ObsType);
[arrQM2, FinalPRNs2, FinalTTs2] = ReadQM(FileQM2);
QM2 = SelectQM(arrQM2, ObsType);

%% 두 QM 파일에서 공통시각(epoch) 추출
FinalTTs = intersect(QM1(:, 1), QM2(:, 1));

%% 항법메시지를 읽어들여서 행렬로 저장하고, Klobuchar 모델 추출
eph = ReadEPH(FileNav);

%% 라이넥스 파일에서 대략적인 관측소 좌표를 뽑아냄
AppPos = TruePosRv; 

%% 추정에 필요한 초기치 설정
MaxIter = 4;
EpsStop = 1e-6;
x = AppPos';


%% 추정과정 시작
NoEpochs = length(FinalTTs);
estm = zeros(NoEpochs, 6);  % : c1(gs), c2/3/4(X/Y/Z of Rv), c5(#Sats_obs), c6(#Sats_used)
nEst = 0;
for j = 1:NoEpochs
    
    gs = FinalTTs(j);
    %% 해당 시각 gs의 관측치 선별 및 공통관측 위성 찾기
    indexQM1 = find(QM1(:,1) == gs);
    QM1eBs = QM1(indexQM1,:);
    indexQM2 = find(QM2(:,1) == gs);
    QM1eRv = QM2(indexQM2,:);
    Sats = intersect(QM1eBs(:, 2), QM1eRv(:, 2));
    NoSats = length(Sats);
    %% 기준위성 RS와 다른위성 OS 선별/ SatsEl - c1(gs), c2(prn), c3(el)
    [SatsEl, indxRS] = PickRSel(gs, Sats, eph, TruePosBs);  % : RS 기준위성
    RS = Sats(indxRS);
    %% 기준위성 좌표 먼저 계산 - Bs 기준
    icol = PickEPH(eph, RS ,gs);
    STT = GetSTTbrdc(gs, RS, eph, TruePosBs);
    tc = gs -STT;
    vec_RS = GetSatPosNC(eph, icol, tc);
    vec_RS = RotSatPos(vec_RS, STT);
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
            %% DD 관측치 생성 파트 --- 추후에 for 루프 밖으로 빼야 함 11/8/14
            OS = Sats(kS);
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
            vec_BsRS = vec_RS - TruePosBs;  com_BsRS = norm(vec_BsRS);
            vec_RvRS = vec_RS - x(1:3)';    com_RvRS = norm(vec_RvRS);
            vec_BsOS = vec_OS - TruePosBs;  com_BsOS = norm(vec_BsOS);
            vec_RvOS = vec_OS - x(1:3)';    com_RvOS = norm(vec_RvOS);
            com = (com_BsRS - com_RvRS) - (com_BsOS - com_RvOS);
            y = obs -com;
            %% H 행렬 계산 파트
            H(1,1) = vec_RvRS(1)/com_RvRS - vec_RvOS(1)/com_RvOS;
            H(1,2) = vec_RvRS(2)/com_RvRS - vec_RvOS(2)/com_RvOS;
            H(1,3) = vec_RvRS(3)/com_RvRS - vec_RvOS(3)/com_RvOS;
            
            HTH = HTH + H'*H;
            HTy = HTy + H'*y;
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
end

%% 측위오차 분석 & 그래프 작성
estm = estm(1:nEst, :);
[dXYZ, dNEV] = PosTErrorsJOON(estm(:, 1), TruePosRv, estm(:, 2:5));
% PosErrorsDD

            