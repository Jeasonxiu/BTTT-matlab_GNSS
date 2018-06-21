%
%   modified by Joonseong Gim, aug 7, 2016
%

clear all; close all;

%% QMfile 호출
FileQM = 'Q160524_ubx1';

%% YY DOY 입력
YY = 16;, DOY = 145;

%% load a PRC file
PRCfile = 'JPRT160524.t1';

%% 불변 변수 설정: 빛의 속도, 관측치
CCC = 299792458.;   % CCC = Speed of Light [m/s]
ObsType = 120;      % 사용할 관측치 설정 - 120: C/A = C1

%% 임계고도각 설정
eleCut = 15;

%% 기준좌표 설정
vrsfile = '160524_1_adm.txt';
VRS = load(vrsfile);
VRS(:,1) = VRS(:,1) +17;
TruePos = VRS(1,5:7);

%% rinex file 입력
obsfile = '160524-ubx1.obs';

%% 만들어진 QMfile을 obsevation Rinex 파일 관측 날짜, 시간으로 변경
rename = renameQMfile(obsfile);
[gw, GD] = ydoy2gwgd(YY, DOY); %: GPS WEEK 결정

%% 사용 파일 입력 및 데이터 추출
% --- GLONASS 방송궤도력 파일 : EphGLO, LeapSecond, TauC
FileNavGLO = strcat('brdc',num2str(DOY),'0.',num2str(YY),'g');
EphGlo = ReadEPH_GLO(FileNavGLO);
TauC = ReadTauC(FileNavGLO); %%
% --- GPS 방송궤도력 파일: 전리층-Klobuchar
FileNav = strcat('brdc',num2str(DOY),'0.',num2str(YY),'n');
LeapSec = GetLeapSec(FileNav); %%
%% GPS QM 파일 읽어들여서 행렬로 저장하고, 사용할 관측치 추출
[arrQM, FinalPRNs, FinalTTs] = ReadQM(FileQM);
QMGPS = SelectQM(arrQM, ObsType);       % GPS C1
QMGLO = SelectQM(arrQM, 320);           % GLO C1
FinalQM = [QMGPS; QMGLO];               % GPS/GLO C1
FinalTTs = intersect(unique(FinalQM(:,1)), VRS(:,1));

[GPSPRC, GLOPRC, PRC_sorted] = PRCsort(PRCfile, FinalQM);

%% 항법메시지를 읽어들여서 행렬로 저장하고, Klobuchar 모델 추출
eph = ReadEPH(FileNav);
[al, be] = GetALBE(FileNav);

%% 항법메시지를 읽어들여서 행렬로 저장하고, Klobuchar 모델 추출
eph = ReadEPH(FileNav);
[al, be] = GetALBE(FileNav);
%% 라이넥스 파일에서 대략적인 관측소 좌표를 뽑아내고 위경도로 변환
% AppPos = [-3041241.741 4053944.143 3859873.640];   % : JPspace B point
AppPos = GetAppPos(obsfile);
if AppPos(1) == 0
    AppPos = TruePos;
end

gd = xyz2gd(AppPos); AppLat = gd(1); AppLon = gd(2);

%% 추정에 필요한 초기치 설정
MaxIter = 10;
EpsStop = 1e-5;
ctr = 1; deltat = 1;

%% 추정과정 시작
NoEpochs = length(FinalTTs);
EstPos = zeros(NoEpochs,5);
nEst = 0;
j=1;
estm = zeros(NoEpochs,6);
visiSat = zeros(NoEpochs,2);
GPS0_c = zeros(NoEpochs,32);

% load('DGPSCPtest.mat')
% TruePos = [-3032234.51900000,4068599.11300000,3851377.46900000];

x = [AppPos ctr]; x = x';
x_d = [AppPos ctr]; x_d = x_d';

for j = 1:NoEpochs
    gs = FinalTTs(j);
    indexQM = find(FinalQM(:,1) == gs);         % gs epoch의 QM
    indexPRC = find(PRC_sorted(:,1) == gs);      % gs epoch의 PRC(GPS/GLO)
    QM1 = FinalQM(indexQM,:);           % GPS/GLO C1 1 Epoch
    PRC1 = PRC_sorted(indexPRC,:);
    NoSats = length(QM1(:,1));
    vec_site = x(1:3)';
    vec_site_d = x_d(1:3)';
    GpsQM = find(QM1(:,3) == 120); GpsSats = length(GpsQM);     % GPS C1
    GloQM = find(QM1(:,3) == 320); GloSats = length(GloQM);     % GLO C1
    visiSat(j,1) = gs; visiSat(j,2) = NoSats; visiSat(j,3) = GpsSats; visiSat(j,4) = GloSats;
    sT= mod(gs,86400)/3600;
    sTh = floor(sT); sTm = sT - sTh;
    
    if j == 1 % 시작시간 위성 위치 array
        %         [Sat_ar] = GetSatPosGLO_new(EphGlo,gs,deltat);
        [Sat_ar] = GetSatPosGLO_my(EphGlo,gs,deltat);
        %         fprintf('시작 시간 위치 계산\n');
    elseif (j ~= 1 && sTm == 0) % 정각일 때
        clear Sat_ar
        %         [Sat_ar] = GetSatPosGLO_new(EphGlo,gs,deltat);
        [Sat_ar] = GetSatPosGLO_my(EphGlo,gs,deltat);
        %         fprintf('%d시 정각 갱신\n',sTh);
    elseif (j ~= 1 && mod(sTm,0.5) == 0) % 30분 일 때
        clear Sat_ar
        %         [Sat_ar] = GetSatPosGLO_new(EphGlo,gs,deltat);
        [Sat_ar] = GetSatPosGLO_my(EphGlo,gs,deltat);
        %         fprintf('%d시 30분 갱신\n',sTh);
    end
    
    ZHD = TropGPTh(vec_site, gw, gs);
    for Iter = 1:MaxIter
        HTH = zeros(4,4);                   % PP
        HTH_d = zeros(4,4);                 % DGNSS
        HTy = zeros(4,1);                   % PP
        HTy_d = zeros(4,1);                 % DGNSS
        
        for i = 1:NoSats
            prn = QM1(i,2); type = QM1(i,3); obs = QM1(i,4);
            if type == 120
                prc = PRC1(find(PRC1(:,2) == prn+100),3);               % DGPS PRC
                icol = PickEPH(eph, prn, gs);
                toe = eph(icol, 8); a = eph(icol, 19); b = eph(icol, 20); c = eph(icol, 21); Tgd = eph(icol, 23);
                %----- 신호전달시간 계산
                STT = GetSTTbrdc(gs, prn, eph, x(1:3)');                % GPS PP
                STT_d = GetSTTbrdc(gs, prn, eph, x_d(1:3)');            % DGPS
                tc = gs - STT;                                          % GPS PP
                tc_d = gs - STT_d;                                        % DGPS
                %----- 위성궤도 계산
                vec_sat = GetSatPosNC(eph, icol, tc);                   % GPS PP
                vec_sat_d = GetSatPosNC(eph, icol, tc_d);               % DGPS
                vec_sat = RotSatPos(vec_sat, STT);                      %: GPS PP 지구자전 고려
                vec_sat_d = RotSatPos(vec_sat_d, STT_d);                  %: DGPS 지구자전 고려
                %----- 최종 RHO 벡터 계산
                vec_rho = vec_sat - x(1:3)';                            % GPS PP
                vec_rho_d = vec_sat_d - x_d(1:3)';                      % DGPS
                rho = norm(vec_rho);                                    % GPS
                rho_d = norm(vec_rho_d);                                % DGPS
                [az,el] = xyz2azel(vec_rho, AppLat, AppLon);
                
                if el >= eleCut %15
%                     W = 1;
                    W = MakeW_elpr(el);
                    
                    dRel = GetRelBRDC(eph, icol, tc);                   % GPS PP
                    dRel_d = GetRelBRDC(eph, icol, tc_d);               % DGPS
                    dtSat = a + b*(tc - toe) + c*(tc - toe)^2 - Tgd + dRel;                 % GPS PP
                    dtSat_d = a + b*(tc_d - toe) + c*(tc_d - toe)^2 - Tgd + dRel_d;         % DGPS
                    dIono = ionoKlob(al, be, gs, az, el, vec_site);                 % Klobuchar
                    dTrop_G = ZHD2SHD(gw, gs, vec_site, el, ZHD);                   % GPT model
                    com = rho + x(4)  - CCC * dtSat + dIono + dTrop_G;              % GPS PP(Klobuchar, GPT Model)
                    com_d = rho_d + x(4)  - CCC * dtSat_d - prc;                    % DGPS
                    y = obs - com;                                          % GPS PP
                    y_d = obs - com_d;                                      % DGPS
                    H = [ -vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 1];                  % GPS PP
                    H_d = [ -vec_rho_d(1)/rho_d -vec_rho_d(2)/rho_d -vec_rho_d(3)/rho_d 1];    % DGPS
                    HTH = HTH + H'*W*H;                         % GPS PP
                    HTH_d = HTH_d + H_d'*W*H_d;                 % DGPS
                    HTy = HTy + H'*W*y;                         % GPS PP
                    HTy_d = HTy_d + H_d'*W*y_d;                 % DGPS
                end
            end
        end
        
        xhat = inv(HTH) * HTy;                              % PP
        xhat_d = inv(HTH_d) * HTy_d;                        % DGNSS
        x = x + xhat;                                       % PP
        x_d = x_d + xhat_d;                                 % DGNSS
        
        if norm(xhat) < EpsStop;
            nEst = nEst + 1;
            estm(nEst,1) = gs;                              % PP
            estm_d(nEst,1) = gs;                            % DGNSS
            estm(nEst,2:5) = x(1:4);                        % PP
            estm_d(nEst,2:5) = x_d(1:4);                      % DGNSS
            fprintf('gs: %6.0f     %2.1f    %2.1f \n',gs,norm(TruePos'-x(1:3)),norm(TruePos'-x_d(1:3)));
            break;
        end
    end
end
estm = estm(1:nEst,:);
%% 측위오차 분석 & 그래프 작성
UBLOX(:,1) = VRS(:,1);
UBLOX(:,2:4) = VRS(:,5:7);
[dXYZ, dNEV] = PosTErrors5(estm, UBLOX, visiSat);
[dXYZ_d, dNEV_d] = PosTErrors5(estm_d, UBLOX, visiSat);
