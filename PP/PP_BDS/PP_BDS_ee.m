clear all; close all

%% test data
% YY = 17; DOY =25; % QMfile = 'QM170125_A'; % navfile = 'brdm0250.17p';
% YY = 17; DOY =25; % QMfile = 'QM170125_B'; % navfile = 'brdm0250.17p';
YY = 15; DOY =11;QMfile = 'QJFNG_15010'; navfile = 'brdm0100.15p';

% TruePos = [-3041235.578 4053941.677 3859881.013];       % A point
% TruePos = [-3041241.741 4053944.143 3859873.640];       % B point
TruePos = [-2279828.8529 5004706.5404 3219777.4631]; %: JFNG site log 121029



%% 불변 변수 설정: 빛의 속도, 관측치
CCC = 299792458.;   % CCC = Speed of Light [m/s]
%% 임계고도각 설정
eleCut = 15;
%% 윤초 핸들링
LeapSecBDS = 14;

%% 사용할 관측치 선정
c_ObsType = 220; % bds C1
c_ObsType_snr = 241;

%% QM 파일 읽어들여서 행렬로 저장하고, 사용할 관측치 추출
[arrQM, FinalPRNs, FinalTTs] = ReadQM(QMfile);

QM = SelectQM(arrQM,c_ObsType);
QM(:,2) = QM(:,2) + 200;        % GPS에서 사용하던 SelectQM 함수 사용시 항법시스템별 구분이 되지 않아 임의로 PRN + 200을 해줌
% QM = QM(find(QM(:,2) ~= 214),:);
FinalTTs = unique(QM(:,1));

%% Ephemeris 행렬 생성
[eph, trashPRN, trashT]=ReadEPH_all(navfile);

%% 항법메시지를 읽어들여서 행렬로 저장하고, Klobuchar 모델 추출
gps_nav = strcat('brdm', num2str(DOY,'%03d'), '0.', num2str(YY,'%02d'), 'p');   %: Navigation RINEX file
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
ctr = 1; 
x = [AppPos ctr]; x = x';

%% 추정과정 시작
NoEpochs = length(FinalTTs);
nEst = 0;

for j = 1:NoEpochs
tic;
    gs = FinalTTs(j);
    indexQM = find(QM(:,1) == gs);
    QM1e = QM(indexQM,:);
    NoSats = length(QM1e);
    
    for iter = 1:Maxiter
        HTH = zeros(4,4);
        HTy = zeros(4,1);
        
        if NoSats <= 4
            break
        end
        vec_site = x(1:3)';
        ZHD = TropGPTh(TruePos, gw, gs); %: TROP: GPT
        NoUsed_bds = 0;
        for i = 1:NoSats
            
            prn = QM1e(i,2);
            obs = QM1e(i,4);
            icol = PickEPH(eph, prn, gs);
            toe = eph(icol, 8); a = eph(icol, 19); b = eph(icol, 20); c = eph(icol, 21); Tgd = eph(icol, 23);
            STT = GetSTTbrdm(gs, eph, icol, x(1:3)); % 신호전달시간 계산
            tc = gs - STT- LeapSecBDS;  % 신호전달 시간 보정... bds윤초 반영
            SatPos = GetSatPosNC_GC(eph, icol, tc); % 위성위치 산출
            SatPos = RotSatPos(SatPos, STT); % 지구자전효과 고려
            vec_rho = SatPos - vec_site';
            rho = norm(vec_rho);
            [az,el] = xyz2azel(vec_rho, AppLat, AppLon);
            
%             dIono = ionoKlob(al, be, gs, az, el, x(1:3)); % 이온층 보정(Klobuchar 모델)
            dIono = klobuchar(al, be, gs, az, el, x(1:3)); % 종석씨가 수정한 이온층 보정(Klobuchar 모델)
            dTrop = ZHD2SHD(gw, gs, TruePos, el, ZHD); % 대류권 보정
            dRel = GetRelBRDC(eph, icol, tc); % 상대성효과
            dtSat = a + b*(tc - toe) + c*(tc - toe)^2 - Tgd + dRel; % 위성시계오차 계산... 그룹딜레이, 상대성효과 보정
            
            if el >=eleCut % 임계고도각
                W = 1;
%                 com = rho + x(4) - CCC*dtSat + dTrop + dIono; % bds 계산값
                com = rho + x(4) - CCC*dtSat + dTrop ; % bds 계산값
                H = [-vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 1];
                y = obs - com;
                HTH = HTH + H'*W*H;
                HTy = HTy + H'*W*y;
                NoUsed_bds = NoUsed_bds + 1;    % 사용 위성 수
            end
        end
        P = inv(HTH);
        xhat = inv(HTH) * HTy;
        x = x + xhat;
%         if norm(xhat) < EpsStop;
            if iter == Maxiter
            nEst = nEst + 1;
            estm(nEst,1) = gs;
            estm(nEst,2:5) = x(1:4);
            estm(nEst,6:7) = [NoSats, NoUsed_bds];
            fprintf('%8d : %3d : %8.2f : %8.2f :%8.2f\n', j, iter, x(1)' - TruePos(1), x(2)' - TruePos(2),  x(3)' - TruePos(3));
            break;
        end
    end
end
toc;
estm = estm(1:nEst,:);
[dXYZ, dNEV] = PosTErrors2(estm(:,1), TruePos, estm(:,2:5),estm(:,6:7));
