%
%   modified by Joonseong Gim, aug 7, 2016
%

clear all; close all;

%% QMfile 호출
FileQM = 'QDRUA_ob153_1';

%% YY DOY 입력
YY = 16;, DOY = 153;

%% load a PRC file
PRCfile = 'JPRT160601.t1';

%% 불변 변수 설정: 빛의 속도, 관측치
CCC = 299792458.;   % CCC = Speed of Light [m/s]
ObsType = 120;      % 사용할 관측치 설정 - 120: C/A = C1

%% 임계고도각 설정
eleCut = 15;

%% 기준좌표 설정
TruePos = [-3041235.578 4053941.677 3859881.013];   % : JPspace A point
% TruePos = [-3041241.741 4053944.143 3859873.640];   % : JPspace B point

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
FinalTTs = unique(FinalQM(:,1));

[GPSPRC, GLOPRC, PRC_sorted] = PRCsort(PRCfile, FinalQM);

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
GLOo_c = zeros(NoEpochs,24);

% load('DGPSCPtest.mat')
% TruePos = [-3032234.51900000,4068599.11300000,3851377.46900000];

x = [AppPos ctr]; x = x';
x_cp = [AppPos ctr]; x_cp = x_cp';


% for j = 1:NoEpochs
for j = 1:908
    gs = FinalTTs(j);
    indexQM = find(FinalQM(:,1) == gs);         % gs epoch의 QM
    indexPRC = find(PRC_sorted(:,1) == gs);      % gs epoch의 PRC(GPS/GLO)
    QM1 = FinalQM(indexQM,:);           % GPS/GLO C1 1 Epoch
    PRC1 = PRC_sorted(indexPRC,:);
    NoSats = length(QM1(:,1));
    vec_site = x(1:3)';
    vec_site_cp = x_cp(1:3)';
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
        HTy = zeros(4,1);                   % PP
        HTy_cp = zeros(4,1);                 % DGNSS
        
        for i = 1:NoSats
            prn = QM1(i,2); type = QM1(i,3); obs = QM1(i,4);
            if type == 120
%                 prc = PRC1(find(PRC1(:,2) == prn+100),3);               % DGPS PRC
%                 if ~isempty(prc)
%                     icol = PickEPH(eph, prn, gs);
%                     toe = eph(icol, 8); a = eph(icol, 19); b = eph(icol, 20); c = eph(icol, 21); Tgd = eph(icol, 23);
%                     %----- 신호전달시간 계산
%                     STT = GetSTTbrdc(gs, prn, eph, x(1:3)');                % GPS PP
%                     tc = gs - STT;                                          % GPS PP
%                     %----- 위성궤도 계산
%                     vec_sat = GetSatPosNC(eph, icol, tc);                   % GPS PP
%                     vec_sat = RotSatPos(vec_sat, STT);                      %: GPS PP 지구자전 고려
%                     %----- 최종 RHO 벡터 계산
%                     vec_rho = vec_sat - x(1:3)';                            % GPS PP
%                     rho = norm(vec_rho);                                    % GPS
%                     [az,el] = xyz2azel(vec_rho, AppLat, AppLon);
%                     
%                     if el >= eleCut %15
%                         W = 1;
% %                         W = MakeW_elpr(el);
%                         
%                         dRel = GetRelBRDC(eph, icol, tc);                   % GPS PP
%                         dtSat = a + b*(tc - toe) + c*(tc - toe)^2 - Tgd + dRel;                 % GPS PP
%                         dIono = ionoKlob(al, be, gs, az, el, vec_site);                 % Klobuchar
%                         dTrop_G = ZHD2SHD(gw, gs, vec_site, el, ZHD);                   % GPT model
%                         com = rho + x(4)  - CCC * dtSat + dIono + dTrop_G;              % GPS PP(Klobuchar, GPT Model)
%                         cp = - dIono - dTrop_G -prc;                            % DGNSS_CP
%                         
%                         y = obs - com;                                          % GPS PP
%                         
%                         H = [ -vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 1 0];                  % GPS PP
%                         HTH = HTH + H'*W*sys_w(1)*H;                         % GPS PP
%                         HTy = HTy + H'*W*sys_w(1)*y;                         % GPS PP
%                         HTy_cp = HTy_cp + H'*W*sys_w(1)*cp;                  % DGPS_CP
%                     end
%                 end
            else
                
                % --- 연직방향 대류권 지연략 계산: 에포크별 한 번 계산
                prc = PRC1(find(PRC1(:,2) == prn+300),3);
                if ~isempty(prc)
                    tc = gs - LeapSec;
                    icol=PickEPH_GLO2(EphGlo, prn, tc);
                    
                    TauN=EphGlo(icol,12); GammaN=EphGlo(icol,13); %: tau & gamma 시계오차 보정에 사용
                    ch_num=EphGlo(icol,16); %: channel number 전리층 보정에 사용
                    
                    % 신호전달시간 계산
                    STT = GetSTTbrdcGLO2(Sat_ar,gs,prn,x(1:3));                 % GLO PP
                    % LeapSecond & 신호전달 시간을 보정한 위성 위치 산출
                    [SatPos, SatVel] = SatPosLS_STT(Sat_ar,gs,prn,LeapSec,STT,TauC);            % GLO PP
                                        
                    % 지구 자전효과 고려
                    SatPos = RotSatPos(SatPos,STT);                     % GLO PP
                    %             SatPos = SatPos';
                    
                    DistXYZ = SatPos - x(1:3)';                         % GLO PP
                    DistNorm = norm(DistXYZ);                           % GLO PP
                    [az,el] = xyz2azel(DistXYZ, AppLat, AppLon);
                    %%
                    if el>=eleCut
                        % 전리층 보정 % 대류권 보정
                        W = 1;
%                         W = MakeW_elpr(el);
                        
                        ttc = tc - TauC;                                    % GLO PP
                        dIono = Klo_R(vec_site,al,be,ttc,SatPos,ch_num);
                        dTrop = ZHD2SHD(gw,gs,vec_site,el,ZHD);
                        % 상대적 효과 (→위성속도 이용)
                        dRel = (-2/CCC^2) * dot(SatPos, SatVel);            % GLO PP
                        % DCB 고려
                        %                     dDCB = AppDCB_glo(DCB,prn-200);
                        %%
                        % 위성시계오차 tsv, tb
                        tsv = tc;                                           % GLO PP
                        tb = EphGlo(icol,2) + LeapSec; % GPStime - 16; / GLOtime + 16; 확인하기
                        %                 dtSat = TauN - GammaN*(tsv-tb) + TauC + dRel + dDCB;
                        dtSat = TauN - GammaN*(tsv-tb) + TauC + dRel;           % GLO PP
                        
                        com = DistNorm + x(4) - CCC*dtSat + dTrop + dIono;      % GLO PP
                        cp = - dTrop - dIono - prc;                             % DGLO_CP 
                        
                        y = obs - com;                                          % GLO PP
                        
                        H = [-DistXYZ(1)/DistNorm -DistXYZ(2)/DistNorm -DistXYZ(3)/DistNorm 1];               % GLO PP
                        HTH = HTH + H'*W*H;                     % GLO PP
                        HTy = HTy + H'*W*y;                     % GLO PP
                        HTy_cp = HTy_cp + H'*W*cp;              % DGLO_cp
                    end
                end
            end
        end
        
        xhat = inv(HTH) * HTy;                              % PP
        xhat_cp = - inv(HTH) * HTy_cp;                        % DGNSS_CP
        x = x + xhat;                                       % PP
        x_cp = x + xhat_cp;                                 % DGNSS_CP
        
        if norm(xhat) < EpsStop;
            nEst = nEst + 1;
            estm(nEst,1) = gs;                              % PP
            estm_cp(nEst,1) = gs;                            % DGNSS_CP
            estm(nEst,2:5) = x(1:4);                        % PP
            estm_cp(nEst,2:5) = x_cp(1:4);                      % DGNSS_CP
            fprintf('gs: %6.0f     %2.1f    %2.1f \n',gs,norm(TruePos'-x(1:3)),norm(TruePos'-x_cp(1:3)));
            break;
        end
    end
end
estm = estm(1:nEst,:);
%% 측위오차 분석 & 그래프 작성
[dXYZ, dNEV] = PosTErrors4(estm(:,1), TruePos, estm(:,2:5),visiSat);
[dXYZ_cp, dNEV_cp] = PosTErrors4(estm_cp(:,1), TruePos, estm_cp(:,2:5),visiSat);
