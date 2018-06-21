function mainGPSkyplotQM(obsfile, AppPos)
%=====================================================
% DO: After reading QM file and then draw skyplot
% Input  : - QMfile=[gs, prn, obstype, measurement]
%          - AppPos=[X, Y, Z]
% Output : figure or jpg file
% by Dong-Hyo Sohn
% (ex) AppPos = [-3026675.978, 4067187.900, 3857246.933];
% (ex) mainGPSkyplotQM('QM15172_ihur', AppPos);
%=====================================================
% clear all; clc; close all;
% QMfile = 'QM15172_ihur';
% AppPos = [-3026675.978, 4067187.900, 3857246.933];

yy = strcat(QMfile(3:4));
doy = strcat(QMfile(5:7));
[yy, doy] = obs2YYDOY(obsfile);
%% �׹��޽��� �о�鿩�� ��ķ� ����
ObsType = 120;
N_file = strcat('brdc',doy,'0.',yy,'n');  % brdc file for GPS
if exist(N_file) == 0
    error('We need a navigation file such as %s',N_file);
    return;
end
eph = ReadEPH(N_file);

%% �뷫���� ������ ��ǥ�� ���浵�� ��ȯ
geod = xyz2gd(AppPos); AppLat = geod(1); AppLon = geod(2);

%% ������-���� ����� �����ϱ� ���� ���ϸ�
% newfile = strcat(QMfile(1:12),'-AzEl.txt');
% fid = fopen(newfile, 'w+');

%% QM ���� �о�鿩�� ��ķ� �����ϰ�, ����� ����ġ ����
[arrQM, FinalPRNs, FinalTTs] = ReadQM(QMfile);
QM = SelectQM(arrQM, ObsType);

%% ������-���� ����ϰ� txt ���Ϸ� �����ϱ�
cnt = 0;
for i = 1: length(FinalTTs)
    if mod(FinalTTs(i),30) ~= 0
        continue;
    end
    
    subQM = [];
    subQM = QM(find(QM(:,1)==FinalTTs(i)),:);
    NoSats = length(subQM);
    gs = subQM(1,1);

    for k = 1 : NoSats
        prn = subQM(k,2);
        icol = PickEPH(eph, prn, gs);           % Pick up the proper column in the ephemerides array
        SatXYZ = GetSatPosNC(eph, icol, gs);    % Compute the XYZ coordinates of a given GPS satellite
        dPos = SatXYZ - AppPos;
        [az,el] = xyz2azel(dPos,AppLat,AppLon); % Compute Azimuth and Elevation
        fprintf(fid,'%8.1f %4d %10.3f %10.3f\n', gs, prn, az, el);
        cnt = cnt + 1;
        QMnew(cnt,:) = [gs, prn, az, el];
    end
end
% fclose(fid);

%% Plot Sky
Az = QMnew(:,3);
El = QMnew(:,4);
cutoff = 5;
site = strcat(QMfile(9:12));

figure(100);
PlotSky(Az, El, cutoff, site);
