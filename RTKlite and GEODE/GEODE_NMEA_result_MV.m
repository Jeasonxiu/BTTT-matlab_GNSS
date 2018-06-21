clear all;
close all;
clc
%% year, month, day
yy = 2016; mo = 12; dd = 07;

%% GEODE & VRS logged file
% geodefile = 'ISMV1321f.txt';
% geodefile = 'ISMV2321g.txt';
% geodefile = 'toDS321g.txt';
% geodefile = 'geode16329.txt';
% geodefile = 'GEODE_test_1.txt';
geodefile = 'GD_16342_T1.txt';
% vrsfile = 'nmea_16321a.txt';
% vrsfile = 'nmea_16321b.txt';
% vrsfile = 'nmea_16321c.txt';
vrsfile = '1207T1.nmea';

%% pick a GPGGA
GPGGA_geode = getGPGGA(geodefile);
filename = geodefile(1,1:max(length(geodefile))-4);
filename = strcat(filename,'_GGA','.txt');
fileoutput = fopen(filename, 'w');

GPGGA_vrs = getGPGGA(vrsfile);
filename_vrs = vrsfile(1,1:max(length(vrsfile))-5);
filename_vrs = strcat(filename_vrs,'_GGA','.txt');
filename_vrs_fix = strcat(filename_vrs,'_GGA_fix','.txt');
filename_vrs_float = strcat(filename_vrs,'_GGA_float','.txt');
filename_vrs_dgps = strcat(filename_vrs,'_GGA_dgps','.txt');
filename_vrs_sa = strcat(filename_vrs,'_GGA_sa','.txt');
fileoutput_vrs = fopen(filename_vrs, 'w');
fileoutput_vrs_fix = fopen(filename_vrs_fix, 'w');
fileoutput_vrs_float = fopen(filename_vrs_float, 'w');
fileoutput_vrs_dgps = fopen(filename_vrs_dgps, 'w');
fileoutput_vrs_sa = fopen(filename_vrs_sa, 'w');

k = 1; ln = 1; fixline = 1; floatline = 1; saline = 1; dgpsline = 1;
for i = 1: length(GPGGA_geode)
    geo_GGA = GPGGA_geode{i,1};
    
    %     line = char(geo_GGA);
    %     fprintf(fileoutput, '%s\r\n', line);
    if length(geo_GGA) > 60
        [hh,mm,ss,x,y,z,la,lo,qi,nSats,ht] = NEWreadGGA3(geo_GGA) ;
        if la > 36 & ht > -10
            [gw, gs] = date2gwgs(yy, mo, dd, hh, mm, ss);
            gpgga_geode(k,:) = [round(gs),x,y,z,la,lo,qi,nSats,ht,ln];
            k = k + 1;
        end
        
    end
    ln = ln + 1;
end
i = 0; k = 1; ln = 1;
for i = 1: length(GPGGA_vrs)
    vrs_GGA = GPGGA_vrs{i,1};
    
    %     line = char(vrs_GGA);
    %     fprintf(fileoutput_vrs, '%s\r\n', line);
    if length(vrs_GGA) > 8
        [hh,mm,ss,x,y,z,la,lo,qi,nSats,ht] = NEWreadGGA3(vrs_GGA) ;
        [gw, gs] = date2gwgs(yy, mo, dd, hh, mm, ss);
        gpgga_vrs(k,:) = [round(gs),x,y,z,la,lo,qi,nSats,ht,ln];
        k = k + 1;
    end
    if qi == 4
        Fixline(fixline,:) = [round(gs),ln];
        fixline = fixline + 1;
    elseif qi == 5
        Floatline(floatline,:) = [round(gs),ln];
        floatline = floatline + 1;
    elseif qi == 1
        SAline(saline,:) = [round(gs),ln];
        saline = saline + 1;
    elseif qi == 2
        DGPSline(dgpsline,:) = [round(gs),ln];
        dgpsline = dgpsline + 1;
    end
    ln = ln + 1;
end

%% VRS, GEODE 공통시간 추출
finalgs = intersect(gpgga_geode(:,1), gpgga_vrs(:,1));

%% 추출된 공통시간에 맞는 NMEA 파일 생성
geode_gga_startend = [max(gpgga_geode(find(gpgga_geode(:,1) == min(finalgs)),10)),...
    max(gpgga_geode(find(gpgga_geode(:,1) == max(finalgs)),10))];

vrs_gga_startend = [max(gpgga_vrs(find(gpgga_vrs(:,1) == min(finalgs)),10)),...
    max(gpgga_vrs(find(gpgga_vrs(:,1) == max(finalgs)),10))];
%% 공통시간 시작, 끝의 gs
Start_gs = gpgga_vrs(vrs_gga_startend(1));
End_gs = gpgga_vrs(vrs_gga_startend(2));
% if min(Fixline(:,1)) > Start_gs
%     vrs_fix_time(1,:) = Fixline(find(Fixline(:,1) == min(Fixline(:,1))),:);
% else
%     vrs_fix_time(1,:) = [Start_gs,min(Fixline(:,2))];
% end
% if max(Fixline(:,1)) > End_gs
%     vrs_fix_time(2,:) = [End_gs,max(Fixline(:,2))];
% else
%     vrs_fix_time(2,:) = Fixline(find(Fixline(:,1) == max(Fixline(:,1))),:);
% end
% Fixline(:,2) = Fixline(:,2) - (vrs_gga_startend(1)-1);
% VRS_GPGGA_fix = GPGGA_vrs(Fixline(:,2)-vrs_gga_startend(1),:);
% VRS_GPGGA_float = GPGGA_vrs(Floatline(:,2)-vrs_gga_startend(1),:);
% VRS_GPGGA_dgps = GPGGA_vrs(DGPSline(:,2)-vrs_gga_startend(1),:);
%% 공통시간에 맞는 GEODE, VRS GPGGA.txt 생성
i = 0; ln = 1;
for i = geode_gga_startend(1):geode_gga_startend(2)
    geo_GGA = GPGGA_geode{i,1};
    line = char(geo_GGA);
    GEODE_GPGGA{ln,1} = line;
    ln = ln + 1;
    fprintf(fileoutput, '%s\r\n', line);
end
i = 0; ln = 1;
for i = vrs_gga_startend(1): vrs_gga_startend(2)
    vrs_GGA = GPGGA_vrs{i,1};
    line = char(vrs_GGA);
    [hh,mm,ss,x,y,z,la,lo,qi,nSats,ht] = NEWreadGGA3(vrs_GGA) ;
    VRS_GPGGA{ln,1} = line;
    if qi == 4
        fprintf(fileoutput_vrs_fix, '%s\r\n', line);
    elseif qi == 5
        fprintf(fileoutput_vrs_float, '%s\r\n', line);
    elseif qi == 1
        fprintf(fileoutput_vrs_sa, '%s\r\n', line);
    elseif qi == 2
        fprintf(fileoutput_vrs_dgps, '%s\r\n', line);
    end
    ln = ln + 1;
    fprintf(fileoutput_vrs, '%s\r\n', line);
end


i = 0;
for i = 1: length(finalgs)
    %     for i = 1: 1
    gs = finalgs(i);
    geo_GGA = gpgga_geode(find(gpgga_geode(:,1) == gs), :);
    vrs_GGA = gpgga_vrs(find(gpgga_vrs(:,1) == gs & gpgga_vrs(:,7) >= 1), :);
    if length(vrs_GGA(:,1)) > 1
        vrs_GGA = vrs_GGA(max(length(vrs_GGA(:,1))),:);
    end
    if ~isempty(vrs_GGA)
        TruePos = vrs_GGA(2:4);
        
        gd = xyz2gd(TruePos); TrueLat = gd(1); TrueLon = gd(2);
        dXYZ = geo_GGA(2:4) - TruePos;
        %% dXYZ를 dNEV로 변환
        dNEV = xyz2topo(dXYZ, TrueLat, TrueLon);
        %% 각 성분별 RMS 계산
        dN = dNEV(:,1); dE = dNEV(:,2); dV = dNEV(:,3);
        dNE = sqrt(dN^2 + dE^2);        %rmsH = myRMS(dNE);
        d2D(i,1) = dNE;
        %rmsV = myRMS(dV);
        d3 = sqrt(dN.^2 + dE.^2 + dV.^2); %rms3 = myRMS(d3);
        d3D(i,1) = d3;
        result(i,:) = [gs, dN, dE, dV, dNE, d3, geo_GGA(7:8), vrs_GGA(7:8)];
    end
end
i = 0; k= 1; kk = 1;
result_dgps = [];
for i = 1: length(finalgs)
    gs = finalgs(i);
    geo_GGA_fix = gpgga_geode(find(gpgga_geode(:,1) == gs), :);
    geo_GGA_dgps = gpgga_geode(find(gpgga_geode(:,1) == gs & gpgga_geode(:,7) == 2) , :);
    vrs_GGA_fix = gpgga_vrs(find(gpgga_vrs(:,1) == gs & gpgga_vrs(:,7) == 4), :);
    if length(vrs_GGA_fix(:,1)) > 1
        vrs_GGA_fix = vrs_GGA_fix(max(length(vrs_GGA_fix(:,1))),:);
    end
    if ~isempty(vrs_GGA_fix)
        TruePos_fix = vrs_GGA_fix(2:4);
        gd = xyz2gd(TruePos_fix); TrueLat = gd(1); TrueLon = gd(2);
        dXYZ = geo_GGA_fix(2:4) - TruePos_fix;
        %% dXYZ를 dNEV로 변환
        dNEV = xyz2topo(dXYZ, TrueLat, TrueLon);
        %% 각 성분별 RMS 계산
        dN = dNEV(:,1); dE = dNEV(:,2); dV = dNEV(:,3);
        dNE = sqrt(dN^2 + dE^2);        %rmsH = myRMS(dNE);
        d2D(k,1) = dNE;
        %rmsV = myRMS(dV);
        d3 = sqrt(dN.^2 + dE.^2 + dV.^2); %rms3 = myRMS(d3);
        d3D(k,1) = d3;
        result_fix(k,:) = [gs, dN, dE, dV, dNE, d3, geo_GGA_fix(7:8), vrs_GGA_fix(7:8)];
        k = k + 1;
    end
    if ~isempty(vrs_GGA_fix) & ~isempty(geo_GGA_dgps)
        TruePos_fix = vrs_GGA_fix(2:4);
        gd = xyz2gd(TruePos_fix); TrueLat = gd(1); TrueLon = gd(2);
        dXYZ = geo_GGA_dgps(2:4) - TruePos_fix;
        %% dXYZ를 dNEV로 변환
        dNEV = xyz2topo(dXYZ, TrueLat, TrueLon);
        %% 각 성분별 RMS 계산
        dN = dNEV(:,1); dE = dNEV(:,2); dV = dNEV(:,3);
        dNE = sqrt(dN^2 + dE^2);        %rmsH = myRMS(dNE);
        d2D(kk,1) = dNE;
        %rmsV = myRMS(dV);
        d3 = sqrt(dN.^2 + dE.^2 + dV.^2); %rms3 = myRMS(d3);
        d3D(kk,1) = d3;
        result_dgps(kk,:) = [gs, dN, dE, dV, dNE, d3, geo_GGA_dgps(7:8), vrs_GGA_fix(7:8)];
        kk = kk + 1;
    end
end

%% rmse & std 계산 all
d2D_all = rms(result(:,5));
d3D_all = rms(result(:,6));
d2D_all_std = std(result(:,5));
d3D_all_std = std(result(:,6));
%% rmse & std 계산 vrs fix
d2D_all_fix = rms(result_fix(:,5));
d3D_all_fix = rms(result_fix(:,6));
d2D_all_std_fix = std(result_fix(:,5));
d3D_all_std_fix = std(result_fix(:,6));
%% rmse & std 계산 vrs fix & geode_dgps
if length(result_dgps) >= 1
    d2D_all_dgps = rms(result_dgps(:,5));
    d3D_all_dgps = rms(result_dgps(:,6));
    d2D_all_std_dgps = std(result_dgps(:,5));
    d3D_all_std_dgps = std(result_dgps(:,6));
end
%% rmse & std display (all)
disp(['d2D all = ', num2str(decimal((d2D_all)*100)/100),'m ', '(',num2str(decimal((d2D_all_std)*100)/100),')'])
disp(['d3D all = ', num2str(decimal((d3D_all)*100)/100),'m ', '(',num2str(decimal((d3D_all_std)*100)/100),')'])
%% rmse & std display (vrs fix)
disp(['d2D all_vrs fix = ', num2str(decimal((d2D_all_fix)*100)/100),'m ', '(',num2str(decimal((d2D_all_std_fix)*100)/100),')'])
disp(['d3D all_vrs fix = ', num2str(decimal((d3D_all_fix)*100)/100),'m ', '(',num2str(decimal((d3D_all_std_fix)*100)/100),')'])
%% rmse & std display (vrs fix)
if length(result_dgps) >= 1
    disp(['d2D all_vrs dgps = ', num2str(decimal((d2D_all_dgps)*100)/100),'m ', '(',num2str(decimal((d2D_all_std_dgps)*100)/100),')'])
    disp(['d3D all_vrs dgps = ', num2str(decimal((d3D_all_dgps)*100)/100),'m ', '(',num2str(decimal((d3D_all_std_dgps)*100)/100),')'])
end

%% all 그래프 plot
tHour = [1:1:length(result(:,1))]';
figure();
subplot(4,4,[1,2,5,6])
plot(result(:,3), result(:,2),'bo'); hold on;
axis([-3 3 -3 3]);
axis square
grid on;
xlabel('\Delta E (meters)')
% xlabel({'\Delta E (meters)',...
%     ['dNE RMSE = ', num2str(decimal(d2D_all)), '   std =', num2str(decimal(d2D_all_std))],...
%     [' 3D RMSE = ', num2str(decimal(d3D_all)), '   std =', num2str(decimal(d3D_all_std))]});
ylabel('\Delta N (meters)')

subplot(4,4,[3,4])
xlim([1 length(tHour)]); grid on; hold on;
plot(tHour(:,1), result(:,2), '.r:', tHour(:,1), result(:,3), '.b:');
legend('\Delta N', '\Delta E')
ylabel('\Delta N,E (meters)');

subplot(4,4,[7,8])
xlim([1 length(tHour)]); grid on; hold on;
plot(tHour(:,1), result(:,4), '.b:'); xlim([1 length(tHour)]); grid on;
legend('\Delta V')
ylabel('\Delta V (meters)')

subplot(4,4,[9,10])
xlim([1 length(tHour)]);
ylim([0 5]);
grid on; hold on;
plot(tHour(:,1), result(:,5), '.b:');
ylabel('Error (meters)')
xlabel('\Delta 2D (meters)')

subplot(4,4,[11,12])
xlim([1 length(tHour)]); grid on; hold on;
plot(tHour(:,1), result(:,7),'.b:');
plot(tHour(:,1), result(:,9),'.r:');
legend('GEODE', 'VRS')
ylabel('Fix Quality');

subplot(4,4,[13,14])
xlim([1 length(tHour)]);
ylim([0 5]);
grid on; hold on;
plot(tHour(:,1), result(:,6), '.b:');
ylabel('Error (meters)')
xlabel('\Delta 3D (meters)')

subplot(4,4,[15,16])
xlim([1 length(tHour)]); grid on; hold on;
plot(tHour(:,1), result(:,8), '.b:');
plot(tHour(:,1), result(:,10), '.r:');
legend('GEODE', 'VRS')
ylabel('Number of Satellites')

%% vrs fix 그래프 plot
tHour_fix = [1:1:length(result_fix(:,1))]';
figure();
subplot(4,4,[1,2,5,6])
plot(result_fix(:,3), result_fix(:,2),'bo'); hold on;
axis([-3 3 -3 3]);
axis square
grid on;
xlabel('\Delta E (meters)')
% xlabel({'\Delta E (meters)',...
%     ['dNE RMSE = ', num2str(decimal(d2D_all)), '   std =', num2str(decimal(d2D_all_std))],...
%     [' 3D RMSE = ', num2str(decimal(d3D_all)), '   std =', num2str(decimal(d3D_all_std))]});
ylabel('\Delta N (meters)')

subplot(4,4,[3,4])
xlim([1 length(tHour_fix)]); grid on; hold on;
plot(tHour_fix(:,1), result_fix(:,2), '.r:', tHour_fix(:,1), result_fix(:,3), '.b:');
legend('\Delta N', '\Delta E')
ylabel('\Delta N,E (meters)');

subplot(4,4,[7,8])
xlim([1 length(tHour_fix)]); grid on; hold on;
plot(tHour_fix(:,1), result_fix(:,4), '.b:'); xlim([1 length(tHour_fix)]); grid on;
legend('\Delta V')
ylabel('\Delta V (meters)')

subplot(4,4,[9,10])
xlim([1 length(tHour_fix)]);
ylim([0 5]);
grid on; hold on;
plot(tHour_fix(:,1), result_fix(:,5), '.b:');
ylabel('Error (meters)')
xlabel('\Delta 2D (meters)')

subplot(4,4,[11,12])
xlim([1 length(tHour_fix)]); grid on; hold on;
plot(tHour_fix(:,1), result_fix(:,7),'.b:');
plot(tHour_fix(:,1), result_fix(:,9),'.r:');
legend('GEODE', 'VRS')
ylabel('Fix Quality');

subplot(4,4,[13,14])
xlim([1 length(tHour_fix)]);
ylim([0 5]);
grid on; hold on;
plot(tHour_fix(:,1), result_fix(:,6), '.b:');
ylabel('Error (meters)')
xlabel('\Delta 3D (meters)')

subplot(4,4,[15,16])
xlim([1 length(tHour_fix)]); grid on; hold on;
plot(tHour_fix(:,1), result_fix(:,8), '.b:');
plot(tHour_fix(:,1), result_fix(:,10), '.r:');
legend('GEODE', 'VRS')
ylabel('Number of Satellites')

if length(result_dgps) >= 1
    [length(result(:,1)), length(result_fix(:,1)), length(result_dgps(:,1))]
else
    [length(result(:,1)), length(result_fix(:,1))]
end

figure()
plot(gpgga_vrs(:,6), gpgga_vrs(:,5), 'r.');
hold on;
plot(gpgga_geode(:,6), gpgga_geode(:,5), 'b.');
plot_google_map;
