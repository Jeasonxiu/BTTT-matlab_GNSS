% function SkyplotQM(QMfile, AppPos)
%=====================================================
% DO: After reading QM file and then draw skyplot
% Input  : - QMfile=[gs, prn, obstype, measurement]
%          - AppPos=[X, Y, Z]
% Output : figure

% Coded by Joonseong GIM, Jan 09, 2017
% (ex) AppPos = [-3026675.978, 4067187.900, 3857246.933];
% (ex) mainGPSkyplotQM('QM15172_ihur', AppPos);
%=====================================================
clear all; clc; close all;
%
%   modified by Joonseong Gim, aug 7, 2016
%

clear all; close all;

%% GPS/GLO 가중치 설정
sys_w = [0.8, 0.2];
%% QMfile 호출
FileQM = 'QDRUA_ob153_1';

%% YY DOY 입력
YY = 16;, DOY = 153;

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
QMGPS = SelectQM(arrQM, 141);       % GPS C1
QMGLO = SelectQM(arrQM, 341);           % GLO C1
FinalQM = [QMGPS; QMGLO];               % GPS/GLO C1
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

%% 추정과정 시작
NoEpochs = length(FinalTTs);
EstPos = zeros(NoEpochs,5);
nEst = 0;
j=1;

% load('DGPSCPtest.mat')
% TruePos = [-3032234.51900000,4068599.11300000,3851377.46900000];

x = [AppPos ctr ctr]; x = x';
x_d = [AppPos ctr ctr]; x_d = x_d';
cntgps = 1; cntglo = 1;

for j = 1:NoEpochs
    % for j = 1:908
    gs = FinalTTs(j);
    sT= mod(gs,86400)/3600;
    sTh = floor(sT); sTm = sT - sTh;
    indexQM = find(FinalQM(:,1) == gs);         % gs epoch의 QM
    QM1 = FinalQM(indexQM,:);           % GPS/GLO C1 1 Epoch
    
    vec_site = x(1:3)';
    GpsQM = QM1(find(QM1(:,3) == 141),:);
    GloQM = QM1(find(QM1(:,3) == 341),:);
    
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
    
    for i = 1:length(GpsQM(:,1))
        prn = GpsQM(i,2); snr = GpsQM(i,4);
        icol = PickEPH(eph, prn, gs);
        vec_sat = GetSatPosNC(eph, icol, gs);                   % GPS PP
        vec_rho = vec_sat - x(1:3)';                            % GPS PP
        [az,el] = xyz2azel(vec_rho, AppLat, AppLon);
        QMnew_gps(cntgps,:) = [gs, prn, az, el, snr];
        cntgps = cntgps + 1;
    end
    
    for k = 1:length(GloQM(:,1))
        prn = GloQM(k,2); snr = GloQM(k,4);
        if~isempty(find(Sat_ar(:,1)==prn & Sat_ar(:,2) == gs))
            tc = gs - LeapSec;
            icol=PickEPH_GLO2(EphGlo, prn, tc);
            % 신호전달시간 계산
            STT = GetSTTbrdcGLO2(Sat_ar,gs,prn,x(1:3));                 % GLO PP
            % LeapSecond & 신호전달 시간을 보정한 위성 위치 산출
            [SatPos, SatVel] = SatPosLS_STT(Sat_ar,gs,prn,LeapSec,STT,TauC);            % GLO PP
            DistXYZ = SatPos - x(1:3)';                         % GLO PP
            [az,el] = xyz2azel(DistXYZ, AppLat, AppLon);
            QMnew_glo(cntglo,:) = [gs, prn, az, el, snr];
            cntglo = cntglo + 1;
        end
    end
end


%% Plot Sky
cutoff = 5;

figure(100);
subplot(3,2,[1,2,3,4])
hold on;
Az = QMnew_gps(:,3); El = QMnew_gps(:,4); Sg = QMnew_gps(:,5);
SgSize = 40;
SgColor = Sg;

%    Draw Circle
circlea = 0:pi/30:2*pi;
cx = cos(circlea);
cy = sin(circlea);

for i= [30 60 90]
    subplot(3,2,[1,2,3,4])
    plot(cx*i, cy*i, '-', 'color', 'k', 'linewidth', 1);
end

for i=[15 45 75]
    subplot(3,2,[1,2,3,4])
    plot(cx*i, cy*i, ':', 'color', 'k', 'linewidth', 1);
    if i == 75
        subplot(3,2,[1,2,3,4])
        plot(cx*i, cy*i, ':', 'color', 'r', 'linewidth', 1);
    end
end

%    Draw Lines inside a Circle
lenmax = 90;
circleTick = (1:6)*pi/6;
cosct = cos(circleTick);
sinct = sin(circleTick);
cax = [-cosct; cosct];
say = [-sinct; sinct];
plot(lenmax*cax, lenmax*say, '-', 'color', 'k', 'linewidth', 1);

xx = (El-90) .* -(sin(Az*pi/180));
yy = (El-90) .* -(cos(Az*pi/180));

%    Draw point Signal Strength value
scatter(xx, yy, SgSize, SgColor,'filled');
title('SkyPlot');
colormap(jet);
colorbar('location', 'EastOutside');
caxis([0 60]);

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

existprn = unique(QMnew_gps(:,2));
for aa = 1:length(existprn)
    prn = existprn(aa);
    if prn < 10
        PRN = strcat('G0',num2str(prn));
    else
        PRN = strcat('G',num2str(prn));
    end
    indxprn = find(QMnew_gps(:,2) == prn);
    QM = QMnew_gps(indxprn,:);
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
Az = QMnew_glo(:,3); El = QMnew_glo(:,4); Sg = QMnew_glo(:,5);
SgSize = 40;
SgColor = Sg;
%    Draw Lines inside a Circle
lenmax = 90;
circleTick = (1:6)*pi/6;
cosct = cos(circleTick);
sinct = sin(circleTick);
cax = [-cosct; cosct];
say = [-sinct; sinct];
plot(lenmax*cax, lenmax*say, '-', 'color', 'k', 'linewidth', 1);

xx = (El-90) .* -(sin(Az*pi/180));
yy = (El-90) .* -(cos(Az*pi/180));

%    Draw point Signal Strength value
scatter(xx, yy, SgSize, SgColor,'filled');
title('SkyPlot');
colormap(jet);
colorbar('location', 'EastOutside');
caxis([0 60]);

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

existprn = unique(QMnew_glo(:,2));
for aa = 1:length(existprn)
    prn = existprn(aa);
    if prn < 10
        PRN = strcat('R0',num2str(prn));
    else
        PRN = strcat('R',num2str(prn));
    end
    indxprn = find(QMnew_glo(:,2) == prn);
    QM = QMnew_glo(indxprn,:);
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

%% SNR bar plot
gpssnrbar = [1:1:32;zeros(1,32)];
glosnrbar = [1:1:27;zeros(1,27)];
% snrbar = zeros(2,59);
for i=1:59
    if i > 32
        if i < 42
            prn = strcat('R','0',num2str(i-32));
            snr = GloQM(find(GloQM(:,2) == i-32),4);
            snrbar{1,i} = prn;
            snrbar{2,i} = snr;
%             SNRbar(1,i) = i;
%             SNRbar(2,i) = snr;
        else
            prn = strcat('R',num2str(i-32));
            snr = GloQM(find(GloQM(:,2) == i-32),4);
            snrbar{1,i} = prn;
            snrbar{2,i} = snr;
%             SNRbar(1,i) = i;
%             SNRbar(2,i) = snr;
        end
    elseif i < 10
        prn = strcat('G','0',num2str(i));
        snr = GpsQM(find(GpsQM(:,2) == i),4);
        snrbar{1,i} = prn;
        snrbar{2,i} = snr;
%         SNRbar(1,i) = i;
%         SNRbar(2,i) = snr;
    else
        prn = strcat('G',num2str(i));
        snr = GpsQM(find(GpsQM(:,2) == i),4);
        snrbar{1,i} = prn;
        snrbar{2,i} = snr;
%         SNRbar(1,i) = i;
%         SNRbar(2,i) = snr;
    end
end
    
for i=1:length(QM1(:,2))
    type=QM1(i,3);
    prn =QM1(i,2);
    snr =QM1(i,4);
    if type < 300
        idxprn = find(gpssnrbar(1,:) == prn);
        gpssnrbar(2,idxprn) = snr;
    else
        idxprn = find(glosnrbar(1,:) == prn);
        glosnrbar(2,idxprn) = snr;
    end
end

figure(100)
subplot(3,2,[5])
bar(gpssnrbar(2,:),'b')
xlim([0 33])
set(gca,'XTick',[1:1:32])
grid on
set(gca,'XTickLabel',{'1', '2','3','4','5','6','7','8','9','10',...
    '11','12','13','14','15','16','17','18','19','20',...
    '21','22','23','24','25','26','27','28','29','30',...
    '31','32'})
xlabel('GPS PRN')
ylabel('SNR');

subplot(3,2,[6])
bar(glosnrbar(2,:),'b')
xlim([0 27])
set(gca,'XTick',[1:1:27])
grid on
set(gca,'XTickLabel',{'1','2','3','4','5','6','7','8','9','10',...
    '11','12','13','14','15','16','17','18','19','20',...
    '21','22','23','24','25','26','27'})
xlabel('GLO PRN')
ylabel('SNR');
