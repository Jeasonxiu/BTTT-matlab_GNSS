%
%   coded by Joonseong Gim, Jan 13, 2016
%

clear all; close all;


%% 불변 변수 설정: 빛의 속도, 관측치
CCC = 299792458.;   % CCC = Speed of Light [m/s]
ObsType = 120;      % 사용할 관측치 설정 - 120: C/A = C1
ObsType2 = 320;      % 사용할 관측치 설정 - 141: snr = S1
TruePos = [-3041235.578 4053941.677 3859881.013];   % : JPspace A point
% TruePos = [-3041241.741 4053944.143 3859873.640];   % : JPspace B point
%% 임계고도각 설정
eleCut = 15;

%% rinex file 입력
% obsfile = 'DDBs115_37.16o';
% navfile = 'DDBs115_37.16n';
% obsfile = '20160211A_2.obs';
% navfile = '20160211A_2.nav';
% obsfile = 'r1.16o';
% navfile = 'r1.16n';
% obsfile = 'l1.obs';
% navfile = 'l1.nav';
% obsfile = 'uatA.obs';
% obsfile = 'SDT1_160824.obs';
obsfile = 'DRUA1533.obs';
% navfile = 'brdc2250.16g';
% obsfile = 'DAUU049r.16.obs';
% navfile = 'brdc0490.16n';
%% Observation Rinex 파일로 부터 QMfile 생성
% WriteObs(obsfile)
% FileQM = 'uatA';
% FileQM = 'QSDT1_16_ob082';
FileQM = 'QDRUA_ob153_3';
%% 만들어진 QMfile을 obsevation Rinex 파일 관측 날짜, 시간으로 변경
rename = renameQMfile(obsfile);
% rename = 'QMDBUU049i_16';
[YY, DOY] = obs2YYDOY(obsfile);
[gw, GD] = ydoy2gwgd(YY, DOY); %: GPS WEEK 결정

%% 사용 파일 입력 및 데이터 추출
% --- GLONASS 방송궤도력 파일 : EphGLO, LeapSecond, TauC
FileNavGLO = strcat('brdc',num2str(DOY),'0.',num2str(YY),'g');
EphGlo = ReadEPH_glo(FileNavGLO);
TauC = ReadTauC(FileNavGLO); %%
% --- GPS 방송궤도력 파일: 전리층-Klobuchar
FileNav = strcat('brdc',num2str(DOY),'0.',num2str(YY),'n');
LeapSec = GetLeapSec(FileNav); %%
%% GPS QM 파일 읽어들여서 행렬로 저장하고, 사용할 관측치 추출
% [arrQM, FinalPRNs, FinalTTs] = ReadQM(rename);
[arrQM, FinalPRNs, FinalTTs] = ReadQM(FileQM);
QMGPS = SelectQM(arrQM, ObsType);       % GPS C1
QMGPS_snr = SelectQM(arrQM, 141);       % GPS snr
QMGLO = SelectQM(arrQM, 320);           % GLO C1
QMGLO_snr = SelectQM(arrQM, 341);       % GLO snr
FinalQM = [QMGPS; QMGLO];               % GPS/GLO C1
FinalQM_snr = [QMGPS_snr; QMGLO_snr];               % GPS/GLO snr
FinalTTs = unique(FinalQM(:,1));
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
x = [AppPos ctr ctr]; x = x';

%% 추정과정 시작
NoEpochs = length(FinalTTs);
EstPos = zeros(NoEpochs,5);
nEst = 0;
j=1;
estm = zeros(NoEpochs,6);
visiSat = zeros(NoEpochs,2);
GPS0_c = zeros(NoEpochs,32);
GLOo_c = zeros(NoEpochs,24);


for j = 1:NoEpochs
    gs = FinalTTs(j);
    indexQM = find(FinalQM(:,1) ==gs);
    QM1 = FinalQM(indexQM,:);           % GPS/GLO C1 1 Epoch
    QM2 = FinalQM_snr(indexQM,:);           % GPS/GLO C1 1 Epoch
    NoSats = length(QM1(:,1));
    vec_site = x(1:3)';
    GpsQM = find(QM1(:,3) == 120); GpsSats = length(GpsQM);     % GPS C1
    GloQM = find(QM1(:,3) == 320); GloSats = length(GloQM);     % GLO C1
    GpsQM_snr = QM2(find(QM2(:,3) == 141),:); GpsSats = length(GpsQM_snr(:,4));     % GPS C1
    GloQM_snr = QM2(find(QM2(:,3) == 341),:); GloSats = length(GloQM_snr(:,4));     % GLO C1
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
        HTH = zeros(5,5);
        HTy = zeros(5,1);
        
        for i = 1:NoSats
            prn = QM1(i,2); type = QM1(i,3); obs = QM1(i,4);
            if type == 120
                S1 = GpsQM_snr(find(GpsQM_snr(:,2) == prn),4);
                icol = PickEPH(eph, prn, gs);
                toe = eph(icol, 8); a = eph(icol, 19); b = eph(icol, 20); c = eph(icol, 21); Tgd = eph(icol, 23);
                %----- 신호전달시간 계산
                STT = GetSTTbrdc(gs, prn, eph, x(1:3)');
                tc = gs - STT;
                %----- 위성궤도 계산
                vec_sat = GetSatPosNC(eph, icol, tc);
                vec_sat = RotSatPos(vec_sat, STT);                      %: 지구자전 고려
                %----- 최종 RHO 벡터 계산
                vec_rho = vec_sat - x(1:3)';
                rho = norm(vec_rho);
                [az,el] = xyz2azel(vec_rho, AppLat, AppLon);
                
                if el >= eleCut %15
                    %                     W = MakeW_elsnr(el,S1);
                    %                     W = MakeW_elpr(el);
                    W = 1;
                    dRel = GetRelBRDC(eph, icol, tc);
                    dtSat = a + b*(tc - toe) + c*(tc - toe)^2 - Tgd + dRel;
                    dIono = ionoKlob(al, be, gs, az, el, vec_site);
                    dTrop_G = ZHD2SHD(gw, gs, vec_site, el, ZHD);                   % GPT model
                    com = rho + x(4)  - CCC * dtSat + dIono + dTrop_G;             % GPT Model
                    y = obs - com;
                    H = [ -vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 1 0];
                    HTH = HTH + H'*W*H;
                    HTy = HTy + H'*W*y;
                end
            else
                
                % --- 연직방향 대류권 지연략 계산: 에포크별 한 번 계산
                S1 = GloQM_snr(find(GloQM_snr(:,2) == prn),4);
                tc = gs - LeapSec;
                icol=PickEPH_GLO2(EphGlo, prn, tc);
                
                TauN=EphGlo(icol,12); GammaN=EphGlo(icol,13); %: tau & gamma 시계오차 보정에 사용
                ch_num=EphGlo(icol,16); %: channel number 전리층 보정에 사용
                
                % 신호전달시간 계산
                STT = GetSTTbrdcGLO2(Sat_ar,gs,prn,x(1:3));
                % LeapSecond & 신호전달 시간을 보정한 위성 위치 산출
                [SatPos, SatVel] = SatPosLS_STT(Sat_ar,gs,prn,LeapSec,STT,TauC);
                
                % 지구 자전효과 고려
                SatPos = RotSatPos(SatPos,STT);
                %             SatPos = SatPos';
                
                DistXYZ = SatPos - x(1:3)';
                DistNorm = norm(DistXYZ);
                [az,el] = xyz2azel(DistXYZ, AppLat, AppLon);
                %%
                if el>=eleCut
                    % 전리층 보정 % 대류권 보정
                    %                     W = MakeW_elsnr(el,S1);
                    %                     W = MakeW_elpr(el);
                    W = 1;
                    ttc = tc - TauC;
                    dIono = Klo_R(vec_site,al,be,ttc,SatPos,ch_num);
                    dTrop = ZHD2SHD(gw,gs,vec_site,el,ZHD);
                    % 상대적 효과 (→위성속도 이용)
                    dRel = (-2/CCC^2) * dot(SatPos, SatVel);
                    % DCB 고려
                    %                     dDCB = AppDCB_glo(DCB,prn-200);
                    %%
                    % 위성시계오차 tsv, tb
                    tsv = tc;
                    tb = EphGlo(icol,2) + LeapSec; % GPStime - 16; / GLOtime + 16; 확인하기
                    %                 dtSat = TauN - GammaN*(tsv-tb) + TauC + dRel + dDCB;
                    dtSat = TauN - GammaN*(tsv-tb) + TauC + dRel;
                    
                    com = DistNorm + x(5) - CCC*dtSat + dTrop + dIono;
                    
                    y = obs - com;
                    H = [-DistXYZ(1)/DistNorm -DistXYZ(2)/DistNorm -DistXYZ(3)/DistNorm 0 1];
                    HTH = HTH + H'*W*H;
                    HTy = HTy + H'*W*y;
                end
            end
        end
        
        xhat = inv(HTH) * HTy;
        x = x + xhat;
        
        if norm(xhat) < EpsStop;
            nEst = nEst + 1;
            estm(nEst,1) = gs;
            estm(nEst,2:5) = x(1:4);
            fprintf('gs: %6.0f     %2.5f \n',gs,norm(TruePos'-x(1:3)));
            break;
        end
    end
end
%         user_gd(j,:) = xyz2gd(estm(j,2:4)); % user의 xyz값 gd 로 변환
%         AppLat = user_gd(j,1); AppLon = user_gd(j,2);
%         user_xyz(j,:) = estm(j,2:4);        % user xyz값 행렬로 변환
%

% estm = estm(find(estm(:,1) > 0),:);
% for es = 1:length(estm(:,1))
%     user_gd(es,:) = xyz2gd(estm(es,2:4)); % user의 xyz값 gd 로 변환
%     AppLat = user_gd(es,1); AppLon = user_gd(es,2);
%     user_xyz(es,:) = estm(es,2:4);        % user xyz값 행렬로 변환
% end
%% user_position mean값 기준 topology 계산 과정
% [user_mean] = PlotMeanTopo(user_xyz);
% user_mean = [mean(user_xyz(:,1)) mean(user_xyz(:,2)) mean(user_xyz(:,3))];
% user_mean = [-3055625.93798130,4034964.65967131,3868139.09930058];
% user_mean = [-3055530.33388512,4035091.49810586,3868085.81846658];
% user_mean = TruePos;
% [dXYZ, dNEV] = PosTErrorsJOON(estm(:, 1), user_mean, estm(:, 2:5));

% [dXYZ, dNEV] = PosTErrorsJOON(estm(:, 1), TruePos, estm(:, 2:5));
% [dXYZ, dNEV] = PosTErrors(estm(:, 1), TruePos, estm(:, 2:5));

%% QM type 별 plot
% % PPPlotQM(rename,141)


%% 구글 plot
% figure(201)
% axis([127.03 127.041 37.528 37.548]);
% plot_google_map;
% axis equal
% axis([127.03 127.041 37.528 37.548]);
%
% for i = 1:length(user_gd)
%     GD_lon(i,1) = user_gd(i,2);
%     GD_la(i,1) = user_gd(i,1);
%     figure(201)
%     grid on
%     hold on
%     plot(GD_lon, GD_la,'ro','markeredgecolor','y','markerfacecolor','r','markersize',3)
%     drawnow
% end

VRS_text = 'SDT1_VRS_16237.txt';
[target,UBLOX] = gapconv2(VRS_text, 0.43);

%% PostErrors
[dXYZ, dNEV] = PosTErrors5(estm, UBLOX, visiSat);

%% PostErrors
% [dXYZ, dNEV] = PosTErrors4(estm(:,1), TruePos, estm(:,2:5),visiSat);

%% 위성별 SNR
% PPPlotQM(QM2,141)

% figure(201)
% grid on
% hold on
% % axis([(min(user_gd(:,2))-0.0001) (max(user_gd(:,2))+0.0001) (min(user_gd(:,1))-0.0001) (max(user_gd(:,1))+0.0001)]);
% plot(user_gd(:,2), user_gd(:,1),'ro','markeredgecolor','y','markerfacecolor','r','markersize',3)
% plot_google_map;