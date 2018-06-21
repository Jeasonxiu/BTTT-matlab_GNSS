% function Skyplotobs(obsfile)
%=====================================================
% DO: After reading QM file and then draw skyplot
% Input  : - QMfile=[gs, prn, obstype, measurement]
%          - AppPos=[X, Y, Z]
% Output : figure or jpg file
% by Dong-Hyo Sohn
% (ex) AppPos = [-3026675.978, 4067187.900, 3857246.933];
% (ex) mainGPSkyplotQM('QM15172_ihur', AppPos);
%=====================================================
clear all; clc; close all;
QMfile = 'QM15172_ihur';
AppPos = [-3026675.978, 4067187.900, 3857246.933];

obsfile = 'SBBs055_14.16o';

%% Obs Rinex 파일로부터 QMfile 
WriteObs(obsfile)
rename = renameQMfile(obsfile);         % Obs 날짜에 맞는 QMfile로 이름 변경
site = obsfile(2);                      % Logging Site

[YY, DOY] = obs2YYDOY(obsfile);
if DOY < 100
    doy = strcat('0',num2str(DOY));
    yy = num2str(YY);
else
    doy = num2str(DOY);
    yy = num2str(YY);
end

navfile = strcat('brdc',doy,'0.',yy,'n');       % brdc file for GPS
if exist(navfile) == 0
    error('We need a navigation file such as %s',navfile);
    return;
end
eph = ReadEPH(navfile);                 % Obs 날짜에 맞는 Brdc파일 읽기

%% Site에 따른 AppPos 좌표 설정
if site == 'A'
    AppPos = [-3041235.578 4053941.677 3859881.013];   % : JPspace A point
    %% 대략적인 관측소 좌표를 위경도로 변환
    geod = xyz2gd(AppPos); AppLat = geod(1); AppLon = geod(2);
    
elseif site == 'B'
    AppPos = [-3041241.741 4053944.143 3859873.640];   % : JPspace B point
    %% 대략적인 관측소 좌표를 위경도로 변환
    geod = xyz2gd(AppPos); AppLat = geod(1); AppLon = geod(2);
else
    User = PPwC(obsfile,navfile);
    AppPos = User(1,2:4);
    geod = xyz2gd(AppPos); AppLat = geod(1); AppLon = geod(2);
%     AppPos = User(:,2:4);
%     for i = 1: length(AppPos(:,1))
%         geod(i,:) = xyz2gd(AppPos(i,:));
%         AppLat(i,:) = geod(i,1);
%         AppLon(i,:) = geod(i,2);
%     end
end
    

%% QM 파일 읽어들여서 행렬로 저장하고, 사용할 관측치 추출
[arrQM, FinalPRNs, FinalTTs] = ReadQM(rename);
QM = SelectQM(arrQM, 120);
QM1 = SelectQM(arrQM, 141);
QM = [QM, QM1(:,4)];
%% 방위각-고도각 계산하고 txt 파일로 저장하기
cnt = 0;
for i = 1: length(FinalTTs)
    if mod(FinalTTs(i),1) ~= 0
        continue;
    end
    
    subQM = [];
    subQM = QM(find(QM(:,1)==FinalTTs(i)),:);
    NoSats = length(subQM(:,1));
    gs = subQM(1,1);

    for k = 1 : NoSats
        prn = subQM(k,2);
        snr = subQM(k,5);
        icol = PickEPH(eph, prn, gs);           % Pick up the proper column in the ephemerides array
        SatXYZ = GetSatPosNC(eph, icol, gs);    % Compute the XYZ coordinates of a given GPS satellite
        dPos = SatXYZ - AppPos;
        [az,el] = xyz2azel(dPos,AppLat,AppLon); % Compute Azimuth and Elevation
%         fprintf(fid,'%8.1f %4d %10.3f %10.3f\n', gs, prn, az, el);
        cnt = cnt + 1;
        QMnew(cnt,:) = [gs, prn, az, el, snr];
    end
end
% fclose(fid);

%% Plot Sky
% Az = QMnew(:,3);
% El = QMnew(:,4);
cutoff = 5;


figure(100);
hold on;
Az = QMnew(:,3); El = QMnew(:,4); Sg = QMnew(:,5);
SgSize = 40;
SgColor = Sg;


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
    
existprn = unique(QMnew(:,2));
for aa = 1:length(existprn)
    prn = existprn(aa);
    if prn < 10
        PRN = strcat('G0',num2str(prn));
    else
        PRN = strcat('G',num2str(prn));
    end
    indxprn = find(QMnew(:,2) == prn);
    QM = QMnew(indxprn,:);
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

% figure(100);
% hold on;
% yy = (El-90).* -(cos(Az*pi/180));
% xx = (El-90).* -(sin(Az*pi/180));
% 
% % 나타내고자하는 방위각, 고도각의 위치에 빨간색 점을 찍는다.
% plot(xx, yy, '.','color', 'r', 'Markersize', 15);
