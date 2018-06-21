function jprtSkyplot(jprtfile)
%=====================================================
% DO: After reading QM file and then draw skyplot
% Input  : jprt prc file
% Output : figure or jpg file

% close all; clear all;

%% prc load
PRC = load(jprtfile);
GPSidx = find(PRC(:,2) < 300);
PRC = PRC(GPSidx,:);
PRC(:,2). - 100;
%% jprt 파일로 부터 YY, doy 생성
yy = str2num(jprtfile(5:6));
mm = str2num(jprtfile(7:8));
dd = str2num(jprtfile(9:10));
doy = date2doy(dd,mm,yy);

%% navfile load
if doy < 100
    navfile = strcat('brdc0',num2str(doy),'0','.',num2str(yy),'n');
else
    navfile = strcat('brdc',num2str(doy),'0','.',num2str(yy),'n');
end
eph = ReadEPH(navfile);

%% NMEAQM 핸들링
FinalTTs = unique(PRC(:,1));

AppPos = [-3026795.499 4067267.161 3857084.459];            %% JPRT 관측소 위치
geod = xyz2gd(AppPos); AppLat = geod(1); AppLon = geod(2);

cnt = 0;
PRN = unique(eph(:,18));
for i = 1: length(FinalTTs)
    if mod(FinalTTs(i),1) ~= 0
        continue;
    end
    
    subQM = [];
    subQM = PRC_Sorted(find(PRC_Sorted(:,1)==FinalTTs(i)),:);
    NoSats = length(subQM);             % PRC 수신 위성
%     NoSats = length(unique(eph(:,18)));             % brdc 전체 위성
    gs = subQM(1,1);

    for k = 1 : NoSats
%         prn = PRN(k,1);                  % brdc 전체 위성
        prn = subQM(k,2);                 % PRC 수신 위성 
        icol = PickEPH(eph, prn, gs);           % Pick up the proper column in the ephemerides array
        SatXYZ = GetSatPosNC(eph, icol, gs);    % Compute the XYZ coordinates of a given GPS satellite
        dPos = SatXYZ - AppPos;
        [az,el] = xyz2azel(dPos,AppLat,AppLon); % Compute Azimuth and Elevation
%         fprintf(fid,'%8.1f %4d %10.3f %10.3f\n', gs, prn, az, el);
        cnt = cnt + 1;
        QMnew(cnt,:) = [gs, prn, az, el];
    end
end
site = ' ';
cutoff = 5;


figure(1000);
Skymap(cutoff, site);

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
    figure(1000);
    hold on;
    plot(xx, yy, '.', 'Markersize', 25);
    plot(xlast, ylast, 'r.', 'Markersize', 25);
    text(xlast, ylast, PRN,'Fontsize',15)
end