%
%   modified by Joonseong Gim, aug 7, 2016
%

clear all; close all;

%% QMfile 호출
FileQM = 'QDRUA_ob153_1';

%% Almanac file & almanac matrix
almfile = '2016_875.txt' ;
[alm] = readAlm_gps(almfile);

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
obsfile = 'DRUA1531.obs';

%% 만들어진 QMfile을 obsevation Rinex 파일 관측 날짜, 시간으로 변경
% rename = renameQMfile_(obsfile);
[gw, GD] = ydoy2gwgd(YY, DOY); %: GPS WEEK 결정

%% 사용 파일 입력 및 데이터 추출
% --- GPS 방송궤도력 파일: 전리층-Klobuchar
FileNav = strcat('brdc',num2str(DOY),'0.',num2str(YY),'n');
LeapSec = GetLeapSec(FileNav); %%
%% GPS QM 파일 읽어들여서 행렬로 저장하고, 사용할 관측치 추출
[arrQM, FinalPRNs, FinalTTs] = ReadQM(FileQM);
QMGPS = SelectQM(arrQM, ObsType);       % GPS C1
QMGPSsnr = SelectQM(arrQM, 141);       % GPS C1
FinalQM = [QMGPS];               % GPS/GLO C1
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

% load('DGPSCPtest.mat')
% TruePos = [-3032234.51900000,4068599.11300000,3851377.46900000];

x = [AppPos ctr]; x = x';
x_alm = [AppPos ctr]; x_alm = x_alm';
cnt = 1;

for j = 1:NoEpochs
    % for j = 1:1
    gs = FinalTTs(j);
    indexQM = find(FinalQM(:,1) == gs);         % gs epoch의 QM
    indexPRC = find(PRC_sorted(:,1) == gs);     % gs epoch의 PRC(GPS)
    indexSNR = find(QMGPSsnr(:,1) == gs);       % gs epoch의 SNR
    QM1 = FinalQM(indexQM,:);                   % GPS C1 1 Epoch
    QM1snr = QMGPSsnr(indexSNR,:);
    PRC1 = PRC_sorted(indexPRC,:);
    NoSats = length(QM1(:,1));
    vec_site = x(1:3)';
    vec_site_alm = x_alm(1:3)';
    GpsQM = find(QM1(:,3) == 120); GpsSats = length(GpsQM);     % GPS C1
    visiSat(j,1) = gs; visiSat(j,2) = NoSats; visiSat(j,3) = GpsSats;
    sT= mod(gs,86400)/3600;
    sTh = floor(sT); sTm = sT - sTh;
    
    ZHD = TropGPTh(vec_site, gw, gs);
%     for Iter = 1:MaxIter
    for Iter = 1:5
        HTH = zeros(4,4);                    % GPS PP
        HTH_alm = zeros(4,4);                % GPS alm PP
        HTy = zeros(4,1);                    % GPS PP
        HTy_alm = zeros(4,1);                % GPS alm PP
        
        for i = 1:NoSats
            prn = QM1(i,2); type = QM1(i,3); obs = QM1(i,4);
            snr = QM1snr(find(QM1snr(:,2)==prn),4);
            
            icol = PickEPH(eph, prn, gs);
            toe = eph(icol, 8); a = eph(icol, 19); b = eph(icol, 20); c = eph(icol, 21); Tgd = eph(icol, 23);
            toa = alm(find(alm(:,1) == prn),4); Af0 = alm(find(alm(:,1) == prn),11); Af1 = alm(find(alm(:,1) == prn),12);
            %----- 신호전달시간 계산
            STT = GetSTTbrdc(gs, prn, eph, x(1:3)');                % GPS PP
            STT_alm = GetSTTalm(gs, prn, alm, x_alm(1:3)');
            tc = gs - STT;                                          % GPS PP
            tc_alm = gs - STT_alm;                                  % GPS alm PP
            %----- 위성궤도 계산
            vec_sat = GetSatPosNC(eph, icol, tc);                   % GPS PP
            vec_sat_alm = GetSatPos_almanac(alm, prn, tc_alm);      % GPS alm PP
            vec_sat = RotSatPos(vec_sat, STT);                      % GPS PP 지구자전 고려
            vec_sat_alm = RotSatPos(vec_sat_alm, STT_alm);          % GPS alm PP 지구자전 고려
            
            %----- 최종 RHO 벡터 계산
            vec_rho = vec_sat' - x(1:3)';                            % GPS PP
            vec_rho_alm = vec_sat_alm' - x_alm(1:3)';                    % GPS alm PP
            
            rho = norm(vec_rho);                                    % GPS
            rho_alm = norm(vec_rho_alm);                            % GPS alm PP
            
            [az,el] = xyz2azel(vec_rho, AppLat, AppLon);
            [az_alm,el_alm] = xyz2azel(vec_rho_alm, AppLat, AppLon);
            
            
            if el >= eleCut %15
                W = 1;
                %                         W = MakeW_elpr(el);
                %                         W = MakeW_elsnr(el,snr);
                
                dRel = GetRelBRDC(eph, icol, tc);                               % GPS PP
                dRel_alm = GetRelAlm(alm, prn, tc_alm);                         % GPS alm PP
                dtSat = a + b*(tc - toe) + c*(tc - toe)^2 - Tgd + dRel;         % GPS PP
                dtSat_alm = Af0 + Af1*(tc_alm - toa) + dRel_alm;                % GPS alm PP
                dIono = ionoKlob(al, be, gs, az, el, vec_site);                 % Klobuchar
                dTrop_G = ZHD2SHD(gw, gs, vec_site, el, ZHD);                   % GPT model
                dIono_alm = ionoKlob(al, be, gs, az_alm, el_alm, vec_site_alm);                 % Klobuchar
                dTrop_G_alm = ZHD2SHD(gw, gs, vec_site_alm, el_alm, ZHD);                   % GPT model
                com = rho + x(4)  - CCC * dtSat + dIono + dTrop_G;              % GPS PP(Klobuchar, GPT Model)
                com_alm = rho_alm + x_alm(4)  - CCC * dtSat_alm + dIono_alm + dTrop_G_alm;          % GPS alm PP(Klobuchar, GPT Model)
                
                y = obs - com;                                          % GPS PP
                y_alm = obs - com_alm;                                  % GPS alm PP
                
                H = [ -vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 1];       % GPS PP
                H_alm = [ -vec_rho_alm(1)/rho_alm -vec_rho_alm(2)/rho_alm...
                    -vec_rho_alm(3)/rho_alm 1];               % DGPS CP alm
                HTH = HTH + H'*W*H;                         % GPS PP
                HTH_alm = HTH_alm + H_alm'*W*H_alm;             % DGPS CP alm
                HTy = HTy + H'*W*y;                         % GPS PP
                HTy_alm = HTy_alm + H_alm'*W*y_alm;                  % DGPS CP alm
                
                if Iter == 1
                    %-----Az, El 저장
                    BRDCazel(cnt,1:4) = [ gs, prn, az,el];
                    Almazel(cnt,1:4) = [gs, prn, az_alm,el_alm];
                    %----- brdc, almanac 위성 좌표 저장
                    SatPos(cnt,1:8) = [gs, prn, vec_sat',vec_sat_alm'];
                    %----- rho 저장
                    Rho(cnt,1:4) = [gs, prn, rho, rho_alm];
                    %----- vec_rho 저장
                    Vec_Rho(cnt,1:8) = [gs, prn, vec_rho, vec_rho_alm];
                    %----- O-C 저장
                    Yvalue(cnt,1:4) = [gs, prn, y, y_alm];
                    %-----H matrix 저장
                    H_matrix(cnt,1:10) = [gs prn H, H_alm];
                    %-----HTy matrix 저장
                    HTy_matrix(cnt,1:10) = [gs, prn, HTy', HTy_alm'];
                    cnt = cnt + 1;
                end
            end
            
        end
        
        xhat = inv(HTH) * HTy;                              % GPS PP
        xhat_alm = inv(HTH_alm) * HTy_alm;                  % GPS alm PP
        if Iter == 1
            XHAT_matrix(j,1:9) = [gs, xhat', xhat_alm'];
        end
        x = x + xhat;                                       % GPS PP
        x_alm = x_alm + xhat_alm;                               % GPS alm PP
        
%         if norm(xhat_alm) < EpsStop;
            if Iter == 5
            nEst = nEst + 1;
            estm(nEst,1) = gs;                              % PP
            estm_alm(nEst,1) = gs;                       % GPS alm PP
            estm(nEst,2:5) = x(1:4);                        % PP
            estm_alm(nEst,2:5) = x_alm(1:4);          % GPS alm PP
            fprintf('gs: %6.0f     %2.1f    %2.1f    \n',gs,norm(TruePos'-x(1:3)),...
                norm(TruePos'-x_alm(1:3)));
%             break;
        end
    end
    %     %---- BRDC CP, Alm CP 저장
    %     brdc_alm(j,1:9) = [gs,xhat_cp', xhat_alm'];
end

estm = estm(1:nEst,:);
estm_alm = estm_alm(1:nEst,:);
[dXYZ, dNEV] = PosTErrors4(estm(:,1), TruePos, estm(:,2:5),visiSat);
[dXYZ_alm, dNEV_alm] = PosTErrors4(estm_alm(:,1), TruePos, estm_alm(:,2:5),visiSat);

figure(100);
hold on;
Az = BRDCazel(:,3); El = BRDCazel(:,4);
SgSize = 20;


%    Draw Circle
circlea = 0:pi/30:2*pi;
cx = cos(circlea);
cy = sin(circlea);

for i= [30 60 90]
    plot(cx*i, cy*i, '-', 'color', 'k', 'linewidth', 0.5);
end

for i=[15 45 75]
    plot(cx*i, cy*i, ':', 'color', 'k', 'linewidth', 0.5);
    if i == 75
        plot(cx*i, cy*i, ':', 'color', 'r', 'linewidth', 0.5);
    end
end

%    Draw Lines inside a Circle
lenmax = 90;
circleTick = (1:6)*pi/6;
cosct = cos(circleTick);
sinct = sin(circleTick);
cax = [-cosct; cosct];
say = [-sinct; sinct];
plot(lenmax*cax, lenmax*say, '-', 'color', 'k', 'linewidth', 0.5);

xx = (El-90) .* -(sin(Az*pi/180));
yy = (El-90) .* -(cos(Az*pi/180));

%    Draw point Signal Strength value
scatter(xx, yy, SgSize,'r');

title('SkyPlot');

%    Insert direction text
rlen = 1.06 * lenmax;
for i = 1 : length(circleTick)
    ticm1 = int2str(i*30);
    ticm2 = int2str(180+i*30);
    if ticm2 == '360'
        ticm2 =' ';
        %ticm2 ='N';
    end
    text( rlen*sinct(i),  rlen * cosct(i), ticm1, 'horizontalalignment', 'center');
    text(-rlen*sinct(i), -rlen * cosct(i), ticm2, 'horizontalalignment', 'center');
end
set(gca,'FontWeight', 'bold');

axis('equal');
axis('off');

existprn = unique(BRDCazel(:,2));
for aa = 1:length(existprn)
    prn = existprn(aa);
    if prn < 10
        PRN = strcat('G0',num2str(prn));
    else
        PRN = strcat('G',num2str(prn));
    end
    indxprn = find(BRDCazel(:,2) == prn);
    QM = BRDCazel(indxprn,:);
    Az = QM(:,3);     El = QM(:,4);
    yy = (El-90).* -(cos(Az*pi/180));
    xx = (El-90).* -(sin(Az*pi/180));
    ylast = yy(length(yy));
    xlast = xx(length(xx));
    figure(100);
    hold on;
    %     plot(xx, yy, '.', 'Markersize', 15);
    text(xlast, ylast, PRN,'Fontsize',15)
end

figure(100);
hold on;
Az = Almazel(:,3); El = Almazel(:,4); 
SgSize = 20;
%    Draw Lines inside a Circle
lenmax = 90;
circleTick = (1:6)*pi/6;
cosct = cos(circleTick);
sinct = sin(circleTick);
cax = [-cosct; cosct];
say = [-sinct; sinct];


xx = (El-90) .* -(sin(Az*pi/180));
yy = (El-90) .* -(cos(Az*pi/180));

%    Draw point Signal Strength value
scatter(xx, yy, SgSize,'b');
title('SkyPlot');

%    Insert direction text
rlen = 1.06 * lenmax;
for i = 1 : length(circleTick)
    ticm1 = int2str(i*30);
    ticm2 = int2str(180+i*30);
    if ticm2 == '360'
        ticm2 =' ';
        %ticm2 ='N';
    end
    text( rlen*sinct(i),  rlen * cosct(i), ticm1, 'horizontalalignment', 'center');
    text(-rlen*sinct(i), -rlen * cosct(i), ticm2, 'horizontalalignment', 'center');
end
set(gca,'FontWeight', 'bold');

axis('equal');
axis('off');

