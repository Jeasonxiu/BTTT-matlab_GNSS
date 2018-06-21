clear all;
close all;
clc

tic
%% YY, DOY 설정
yy = 2016; mo = 11;, dd = 17;
doy = 322;

%% QMfile load
baseqm = 'rtklite_base_161117';
roverqm = 'rtklite_rover_161117';

%% navfile load
navfile = strcat('brdc',num2str(doy),'0.',num2str(yy-2000),'n');       % brdc file for GPS
if exist(navfile) == 0
    error('We need a navigation file such as %s',navfile);
    return;
end
eph = ReadEPH(navfile);                 % Obs 날짜에 맞는 Brdc파일 읽기

%% base QM 파일 읽어들여서 행렬로 저장하고, 사용할 관측치 추출
[arrQM_base, FinalPRNs_base, FinalTTS_base] = ReadQM(baseqm);
QM_base = SelectQM(arrQM_base, 120);
QM1_base = SelectQM(arrQM_base, 141);
QM_base = [QM_base, QM1_base(:,4)];

%% rover QM 파일 읽어들여서 행렬로 저장하고, 사용할 관측치 추출
[arrQM_rover, FinalPRNs_rover, FinalTTS_rover] = ReadQM(roverqm);
QM_rover = SelectQM(arrQM_rover, 120);
QM1_rover = SelectQM(arrQM_rover, 141);
QM_rover = [QM_rover, QM1_rover(:,4)];

%% RTKLITE logged file
% filename = 'nmea_20161113.txt';
filename = 'rtklite_161117.txt';

% %% 대략좌표
% site = [-3027409.89954179 4047618.47439957 3877061.37645965];                      % Logging Site
% AppPos = [-3027409.89954179 4047618.47439957 3877061.37645965];                      % Logging Site
% geod = xyz2gd(AppPos); AppLat = geod(1); AppLon = geod(2);

%% NMEA text에서 GPGGA, GPGSA, GNGGA 추출
GPGGA = getGPGGA(filename);
GPGSA = getGPGSA(filename);

GN = getGNGGA(filename);

GNGGA_issue = [];
GPGGA_issue = [];
% 
%% NMEA 결과 분석 시작
gp = 1; gn = 1;
usedGPS = zeros(length(GPGGA),13);

for i = 1:length(GPGGA)
    GGA = GPGGA{i,1};
    %% GPGSA가 존재하면 사용된 GPS 위성 수 저장
    if ~isempty(GPGSA)
        GSA = strsplit(GPGSA{i,1},',');
        
        if str2num(cell2mat(GSA(3))) == 3
            num_gps(i,1) = length(GSA) - 7;
            for s = 4:length(GSA) - 4;
                usedGPS(i,s-2) = str2num(cell2mat(GSA(s)));
            end
        else
            num_gps(i,1) = 0;
        end
    else
        num_gps(i,1) = 0;
    end 
    
    if GGA(1:6) == '$GPGGA'
        GNGGA(i,1) = 1;
    elseif GGA(1:6) == '$GNGGA'
        GNGGA(i,1) = 2;
    end
    [hh,mm,ss,x,y,z,la,lo,qi,nSats,ht] = NEWreadGGA3(GGA) ;
    [gw, gs] = date2gwgs(yy, mo, dd, hh, mm, ss);
    GNGGA(i,2:13) = [hh,mm,ss,x,y,z,la,lo,qi,nSats,ht,round(gs)];
    usedGPS(i,1) = round(gs);
end

%% Reference 좌표와 대략좌표 
TruePos = [mean(GNGGA(find(GNGGA(:,10)==4),5)), mean(GNGGA(find(GNGGA(:,10)==4),6)), mean(GNGGA(find(GNGGA(:,10)==4),7))];
AppPos = TruePos;
gd = xyz2gd(TruePos); TrueLat = gd(1); TrueLon = gd(2);
geod = xyz2gd(AppPos); AppLat = geod(1); AppLon = geod(2);

i=0;
for i=1:length(GNGGA)
    GGA = GNGGA(i,:);
    dXYZ = [GGA(5), GGA(6), GGA(7)] - TruePos;
    %% dXYZ를 dNEV로 변환
    dNEV = xyz2topo(dXYZ, TrueLat, TrueLon);
    %% 각 성분별 RMS 계산
    dN = dNEV(:,1); dE = dNEV(:,2); dV = dNEV(:,3);
    dNE = sqrt(dN^2 + dE^2);        %rmsH = myRMS(dNE);
    d2D(i,1) = dNE; d2D(i,2) = GNGGA(i,1); d2D(i,3) = GGA(10);
    %rmsV = myRMS(dV);
    d3 = sqrt(dN.^2 + dE.^2 + dV.^2); %rms3 = myRMS(d3);
    d3D(i,1) = d3; d3D(i,2) = GNGGA(i,1); d3D(i,3) = GGA(10);
    result(i,:) = [dN, dE, dV, dNE, d3, GGA(10), num_gps(i,1), GGA(11), GNGGA(i,1), i, GGA(13)];
%     if GNGGA(i,1) == 1 & qi > 4
%         GPGGA_issue(gp,1) = i;
%         GPGGA_issue(gp,2:10) = result(i,:);
%         gp= gp + 1;
%     elseif GNGGA(i,1) == 2
%         GNGGA_issue(gn,1) = i;
%         GNGGA_issue(gn,2:10) = result(i,:);
%         gn = gn + 1;
%     end
end

%% 공통시간
FinalTTs_base = intersect(FinalTTS_base(:,1), GNGGA(:,13)); 
FinalTTs_rover = intersect(FinalTTS_rover(:,1), GNGGA(:,13)); 
FinalTTs = intersect(FinalTTs_base, FinalTTs_rover);
%% base 방위각-고도각 계산하고 txt 파일로 저장하기
i=0; cnt = 0; cnt2 = 0; cnt3 = 0;
for i = 1: length(FinalTTs)
    if mod(FinalTTs(i),1) ~= 0
        continue;
    end
    
    subQM_base = []; subQM_rover = [];
    subQM_base = QM_base(find(QM_base(:,1)==FinalTTs(i)),:);
    subQM_rover = QM_rover(find(QM_rover(:,1)==FinalTTs(i)),:);
    usedgps = usedGPS(find(usedGPS(:,1)==FinalTTs(i)),:);
    NoSats_base = length(subQM_base);
    NoSats_rover = length(subQM_rover);
    used_GPS = intersect(subQM_rover(:,2), usedgps(:,2:13));
    NoSats_base_used = length(used_GPS);
    gs_base = subQM_base(1,1);

    for k = 1 : NoSats_base
        prn_base = subQM_base(k,2);
        snr_base = subQM_base(k,5);
        icol_base = PickEPH(eph, prn_base, gs_base);           % Pick up the proper column in the ephemerides array
        SatXYZ_base = GetSatPosNC(eph, icol_base, gs_base);    % Compute the XYZ coordinates of a given GPS satellite
        dPos_base = SatXYZ_base - AppPos;
        [az_base,el_base] = xyz2azel(dPos_base,AppLat,AppLon); % Compute Azimuth and Elevation
%         fprintf(fid,'%8.1f %4d %10.3f %10.3f\n', gs, prn, az, el);
        cnt = cnt + 1;
        QMnew_base(cnt,:) = [gs_base, prn_base, az_base, el_base, snr_base];
    end
    for kk = 1 : NoSats_rover
        prn_rover = subQM_rover(kk,2);
        snr_rover = subQM_rover(kk,5);
        icol_rover = PickEPH(eph, prn_rover, gs_base);           % Pick up the proper column in the ephemerides array
        SatXYZ_rover = GetSatPosNC(eph, icol_rover, gs_base);    % Compute the XYZ coordinates of a given GPS satellite
        dPos_rover = SatXYZ_rover - AppPos;
        [az_rover,el_rover] = xyz2azel(dPos_rover,AppLat,AppLon); % Compute Azimuth and Elevation
        %         fprintf(fid,'%8.1f %4d %10.3f %10.3f\n', gs, prn, az, el);
        cnt2 = cnt2 + 1;
        QMnew_rover(cnt2,:) = [gs_base, prn_rover, az_rover, el_rover, snr_rover];
    end
    for kkk = 1:NoSats_base_used
        prn_used = used_GPS(kkk);
        snr_used = subQM_rover(find(subQM_rover(:,2) == prn_used), 5);
        icol_used = PickEPH(eph, prn_used, gs_base);           % Pick up the proper column in the ephemerides array
        SatXYZ_used = GetSatPosNC(eph, icol_used, gs_base);    % Compute the XYZ coordinates of a given GPS satellite
        dPos_used = SatXYZ_used - AppPos;
        [az_used,el_used] = xyz2azel(dPos_used,AppLat,AppLon); % Compute Azimuth and Elevation
%         fprintf(fid,'%8.1f %4d %10.3f %10.3f\n', gs, prn, az, el);
        cnt3 = cnt3 + 1;
        QMnew_used(cnt3,:) = [gs_base, prn_used, az_used, el_used, snr_used];
    end
end
% fclose(fid);
toc
%% Plot Sky
% Az = QMnew(:,3);
% El = QMnew(:,4);
cutoff = 5;

%% base skyplot : figre(100), figure(101)
figure(100);
hold on;
Az_base = QMnew_base(:,3); El_base = QMnew_base(:,4); Sg_base = QMnew_base(:,5);
SgSize_base = 40;
SgColor_base = Sg_base;


%    Draw Circle
    circlea = 0:pi/30:2*pi;
    cx = cos(circlea);
    cy = sin(circlea);
 
    for i= [30 60 90]
        plot(cx*i, cy*i, '-', 'color', 'k', 'linewidth', 1);
    end
 
    for i=[15 45 75]
        plot(cx*i, cy*i, ':', 'color', 'k', 'linewidth', 1);
    end

%    Draw Lines inside a Circle
    lenmax = 90;
    circleTick = (1:6)*pi/6;
    cosct = cos(circleTick); 
    sinct = sin(circleTick);
    cax = [-cosct; cosct];
    say = [-sinct; sinct];
    plot(lenmax*cax, lenmax*say, '-', 'color', 'k', 'linewidth', 1);

    xx_base = (El_base-90) .* -(sin(Az_base*pi/180));
    yy_base = (El_base-90) .* -(cos(Az_base*pi/180));

%    Draw point Signal Strength value
    scatter(xx_base, yy_base, SgSize_base, SgColor_base,'filled');
    title('SkyPlot(base)');
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
    
existprn_base = unique(QMnew_base(:,2));
for aa = 1:length(existprn_base)
    prn_base = existprn_base(aa);
    if prn_base < 10
        PRN_base = strcat('G0',num2str(prn_base));
    else
        PRN_base = strcat('G',num2str(prn_base));
    end
    indxprn_base = find(QMnew_base(:,2) == prn_base);
    QM_base = QMnew_base(indxprn_base,:);
    Az_base = QM_base(:,3);     El_base = QM_base(:,4);
    yy_base = (El_base-90).* -(cos(Az_base*pi/180));
    xx_base = (El_base-90).* -(sin(Az_base*pi/180));
    ylast_base = yy_base(length(yy_base));
    xlast_base = xx_base(length(xx_base));
    figure(100);
    hold on;
%     plot(xx, yy, '.', 'Markersize', 15);
    text(xlast_base, ylast_base, PRN_base,'Fontsize',15)
end


figure(101);
subplot(4,1,[1])
hold on;

hold on;
Az_used = QMnew_used(:,3); El_used = QMnew_used(:,4); Sg_used = QMnew_used(:,5);


SgSize_used = 40;
SgColor_used = Sg_used;
%    Draw point Signal Strength value
% plot(QMnew_used(:,1),15,'r-')
scatter(QMnew_used(:,1), El_used, SgSize_used, SgColor_used,'filled');
xlim([min(QMnew_used(:,1)), max(QMnew_used(:,1))])
title('SkyPlot');
colormap(jet);
% colorbar('location', 'EastOutside');
caxis([0 60]);

aa=0;
existprn_used = unique(QMnew_used(:,2));
for aa = 1:length(existprn_used)
    prn_used = existprn_used(aa);
    if prn_used < 10
        PRN_used = strcat('G0',num2str(prn_used));
    else
        PRN_used = strcat('G',num2str(prn_used));
    end
    indxprn_used = find(QMnew_used(:,2) == prn_used);
    QM_used = QMnew_used(indxprn_used,:);
%     prn_used_last = [max(find(QMnew_used(:,1) == max(QM_used(:,1)))), QM_used(find(QM_used(:,1) == max(QM_used(:,1))),4)];
    prn_used_last = [max(QM_used(:,1)), QM_used(find(QM_used(:,1) == max(QM_used(:,1))),4)];
   
    figure(101);
    subplot(4,1,[1])
    hold on;
    plot([0,max(result(:,11))],[15,15],'r-')
    text(prn_used_last(1), prn_used_last(2), PRN_used,'Fontsize',15);
end

load('rtklite_bds_num');
bdsgs = X(1,:)';
bdsnum = Y(1,:)';

subplot(4,1,[3])
xlim([min(QMnew_used(:,1)), max(QMnew_used(:,1))]); grid on; hold on;
plot(result(find(result(:,6) == 4),11), result(find(result(:,6) == 4),8), 'b.:'); 
plot(result(:,11), result(:,7), 'r.:');
plot(bdsgs(:,1), bdsnum(:,1), 'g.:');
plot(result(find(result(:,6) == 5),11), result(find(result(:,6) == 5),8), 'k*'); 
plot(result(find(result(:,6) == 1),11), result(find(result(:,6) == 1),8), 'r^'); 
% legend('Num of all Sats', 'Num of GPS')
ylabel('Number of Satellites')

subplot(4,1,[4])
xlim([min(QMnew_used(:,1)), max(QMnew_used(:,1))]); grid on; hold on;
plot(result(find(result(:,6) == 4),11), result(find(result(:,6) == 4),4), '.b:');
plot(result(find(result(:,6) == 5),11), result(find(result(:,6) == 5),4), 'k*');
plot(result(find(result(:,6) == 1),11), result(find(result(:,6) == 1),4), 'r^:');
ylim([0 0.05]);
ylabel('Error (meters)')
legend('\Delta NE')
% figure(100);
% hold on;
% yy = (El-90).* -(cos(Az*pi/180));
% xx = (El-90).* -(sin(Az*pi/180));
% 
% % 나타내고자하는 방위각, 고도각의 위치에 빨간색 점을 찍는다.
% plot(xx, yy, '.','color', 'r', 'Markersize', 15);
