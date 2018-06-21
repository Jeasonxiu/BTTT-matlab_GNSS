clear all;
close all;
clc

%% year, month, day
yy = 2016; mo = 12; dd = 09;

%% GEODE logged file
% filename = 'DSB2316i.txt';
% filename = 'DSA2316i.txt';
% filename = 'DSB1316i.txt';
% filename = 'DSA1316i.txt';
% filename = 'DSC1320g.txt';
% filename = 'DSC2320g.txt';
% filename = 'DSD1320h.txt';
% filename = 'GEODE_test_toDS.txt';
filename = '1209_1621_C_GEODE';
% filename = '1209_1554_D_GEODE';
% filename = '1209_1539_E_GEODE';
% filename = '1209_1509_F_GEODE';
% filename = 'DSDC320h.txt';
% filename = 'nmea_20161113.txt';
%% Reference Point
% TruePos = [-3041235.578 4053941.677 3859881.013];   % : JPspace A point
% TruePos = [-3041241.741 4053944.143 3859873.640];   % : JPspace B point
TruePos = [-3041210.419 4053863.515 3859778.262];   % : JPspace C point;
% TruePos = [-3041230.128 4053871.023 3859754.445];   % : JPspace D point
% TruePos = [-3041162.756 4053998.170 3859674.940];   % : JPspace E point
% TruePos = [-3041092.515 4054082.611 3859643.389];   % : JPspace F point
% TruePos = [-3027413.186	4047615.838	3877061.269];

gd = xyz2gd(TruePos); TrueLat = gd(1); TrueLon = gd(2);

GPGGA = getGPGGA(filename);
GN = getGNGGA(filename);

%
dgps = 1; sa = 1;
GPGGA_sa = [];
for i = 1:length(GPGGA)
    GGA = GPGGA{i,1};
    if GGA(1:6) == '$GPGGA'
        GNGGA(i,1) = 1;
    elseif GGA(1:6) == '$GNGGA'
        GNGGA(i,1) = 2;
    end
    if length(GGA) > 60
        [hh,mm,ss,x,y,z,la,lo,qi,nSats,ht] = NEWreadGGA3(GGA) ;
        [gw, gs] = date2gwgs(yy, mo, dd, hh, mm, ss);
        GNGGA(i,1:12) = [round(gs), hh,mm,ss,x,y,z,la,lo,qi,nSats,ht];
        
        dXYZ = [x, y, z] - TruePos;
        %% dXYZ를 dNEV로 변환
        dNEV = xyz2topo(dXYZ, TrueLat, TrueLon);
        %% 각 성분별 RMS 계산
        dN = dNEV(:,1); dE = dNEV(:,2); dV = dNEV(:,3);
        dNE = sqrt(dN^2 + dE^2);        %rmsH = myRMS(dNE);
        d2D(i,1) = dNE; d2D(i,2) = GNGGA(i,1);
        %rmsV = myRMS(dV);
        d3 = sqrt(dN.^2 + dE.^2 + dV.^2); %rms3 = myRMS(d3);
        d3D(i,1) = d3; d3D(i,2) = GNGGA(i,1);
        result(i,:) = [dN, dE, dV, dNE, d3, qi, nSats, GNGGA(i,1), i];
        if qi == 2
            GPGGA_dgps(dgps,1) = i;
            GPGGA_dgps(dgps,2:10) = result(i,:);
            dgps= dgps + 1;
        elseif qi == 1
            GPGGA_sa(sa,1) = i;
            GPGGA_sa(sa,2:10) = result(i,:);
            sa= sa + 1;
            
        end
    end
end
tHour = [1:1:length(result(:,1))]';
if ~isempty(GPGGA_sa)
    d2D_all = rms(d2D(:,1)); d3D_all = rms(d3D(:,1));
    d2D_dgps = rms(GPGGA_dgps(:,5));
    d3D_dgps = rms(GPGGA_dgps(:,6));
    d2D_sa = rms(GPGGA_sa(:,5));
    d3D_sa = rms(GPGGA_sa(:,6));
    d2D_all_std = std(d2D(:,1))
    d3D_all_std = std(d3D(:,1))
    d2D_dgps_std = std(GPGGA_dgps(:,5));
    d3D_dgps_std = std(GPGGA_dgps(:,6));
    d2D_sa_std = std(GPGGA_sa(:,5));
    d3D_sa_std = std(GPGGA_sa(:,5));
    disp(['d2D all = ', num2str(decimal((d2D_all)*100)/100),'m ', '(',num2str(decimal((d2D_all_std)*100)/100),')'])
    disp(['d3D all = ', num2str(decimal((d3D_all)*100)/100),'m ', '(',num2str(decimal((d3D_all_std)*100)/100),')'])
    disp(['d2D DGPS = ', num2str(decimal((d2D_dgps)*100)/100),'m ', '(',num2str(decimal((d2D_dgps_std)*100)/100),')'])
    disp(['d3D DGPS = ', num2str(decimal((d3D_dgps)*100)/100),'m ', '(',num2str(decimal((d3D_dgps_std)*100)/100),')'])
    disp(['d2D SA = ', num2str(decimal((d2D_sa)*100)/100),'m ', '(',num2str(decimal((d2D_sa_std)*100)/100),')'])
    disp(['d3D SA = ', num2str(decimal((d3D_sa)*100)/100),'m ', '(',num2str(decimal((d3D_sa_std)*100)/100),')'])
else
    d2D_all = rms(d2D(:,1)); d3D_all = rms(d3D(:,1));
    %     d2D_under = rms(GPGGA_good(:,5));
    %     d3D_under = rms(GPGGA_good(:,6));
    d2D_dgps = rms(GPGGA_dgps(:,5));
    d3D_dgps = rms(GPGGA_dgps(:,6));
    d2D_all_std = std(d2D(:,1))
    d3D_all_std = std(d3D(:,1))
    %     d2D_under_std = std(GPGGA_good(:,5));
    %     d3D_under_std = std(GPGGA_good(:,6));
    d2D_dgps_std = std(GPGGA_dgps(:,5));
    d3D_dgps_std = std(GPGGA_dgps(:,5));
    % error_result = [decimal((d2D_all)*100)/100, decimal((d3D_all)*100)/100;...
    %     decimal((d2D_GPGGA)*100)/100, decimal((d3D_GPGGA)*100)/100;...
    %     decimal((d2D_GNGGA)*100)/100, decimal((d3D_GNGGA)*100)/100];
    disp(['d2D all = ', num2str(decimal((d2D_all)*100)/100),'m ', '(',num2str(decimal((d2D_all_std)*100)/100),')'])
    disp(['d3D all = ', num2str(decimal((d3D_all)*100)/100),'m ', '(',num2str(decimal((d3D_all_std)*100)/100),')'])
    % disp(['d2D Under = ', num2str(decimal((d2D_under)*100)/100),'m ', '(',num2str(decimal((d2D_under_std)*100)/100),')'])
    % disp(['d3D Under = ', num2str(decimal((d3D_under)*100)/100),'m ', '(',num2str(decimal((d3D_under_std)*100)/100),')'])
    disp(['d2D DGPS = ', num2str(decimal((d2D_dgps)*100)/100),'m ', '(',num2str(decimal((d2D_dgps_std)*100)/100),')'])
    disp(['d3D DGPS = ', num2str(decimal((d3D_dgps)*100)/100),'m ', '(',num2str(decimal((d3D_dgps_std)*100)/100),')'])
end
figure();
hold on;

subplot(4,4,[1,2,5,6])
if ~isempty(GPGGA_sa)
    plot(GPGGA_sa(:,3), GPGGA_sa(:,2),'ro'); hold on;
    plot(GPGGA_dgps(:,3), GPGGA_dgps(:,2),'bo');
    lineCross(0,0,'r',0.1)
    legend('\Delta NE StandAlone', '\Delta NE DGPS(SBAS)')
else
    plot(GPGGA_dgps(:,3), GPGGA_dgps(:,2),'bo');
    lineCross(0,0,'r',0.1)
end
hold on;
% axis([-3 3 -3 3]);
axis square
grid on;
% xlabel({'\Delta E (meters)',...
%     ['dNE RMSE = ', num2str(decimal(rms(result(:,4)))), '   std =', num2str(decimal(std(result(:,4))))],...
%     [' 3D RMSE = ', num2str(decimal(rms(result(:,5)))), '   std =', num2str(decimal(std(result(:,5))))]});
% xlabel({'\Delta E (meters)',...
%     ['dNE RMSE = ', num2str(decimal(rms(result(:,4)))), '   std =', num2str(decimal(std(result(:,4)))),...
%     ' 3D RMSE = ', num2str(decimal(rms(result(:,5)))), '   std =', num2str(decimal(std(result(:,5))))]});
; ylabel('\Delta N (meters)')
%% 그래프 우측
subplot(4,4,[3,4])
plot(tHour(:,1), result(:,1), '.r:', tHour(:,1), result(:,2), '.b:'); xlim([1 length(tHour)]); grid on;
legend('\Delta N', '\Delta E')
ylabel('\Delta N,E (meters)');
subplot(4,4,[7,8])
xlim([1 length(tHour)]); grid on; hold on;
if ~isempty(GPGGA_sa)
    plot(GPGGA_sa(:,1), GPGGA_sa(:,4), '.r');
    plot(GPGGA_dgps(:,1), GPGGA_dgps(:,4), '.b');
    legend('\Delta V(SA)', '\Delta V(DGPS)')
else
    plot(GPGGA_dgps(:,1), GPGGA_dgps(:,4), '.b');
end
ylabel('\Delta V (meters)')

subplot(4,4,[9,10])
xlim([1 length(tHour)]);
ylim([0 5]);
grid on; hold on;
if ~isempty(GPGGA_sa)
    plot(GPGGA_sa(:,1), GPGGA_sa(:,5), '.r');
    plot(GPGGA_dgps(:,1), GPGGA_dgps(:,5), '.b');
    legend('\Delta NE(SA)','\Delta NE(DGPS)')
else
    plot(GPGGA_dgps(:,1), GPGGA_dgps(:,5), '.b');
end
ylabel('2D Error (meters)')

subplot(4,4,[11,12])
plot(tHour(:,1), result(:,6),'.b:'); xlim([1 length(tHour)]); grid on;
ylabel('Fix Quality');
subplot(4,4,[13,14])
xlim([1 length(tHour)]);
ylim([0 5]);
grid on; hold on;
if ~isempty(GPGGA_sa)
    plot(GPGGA_sa(:,1), GPGGA_sa(:,6), '.r');
    plot(GPGGA_dgps(:,1), GPGGA_dgps(:,6), '.b');
    legend('\Delta 3D(SA)','\Delta 3D(DGPS)')
else
    plot(GPGGA_dgps(:,1), GPGGA_dgps(:,6), '.b');
end
ylabel('3D Error (meters)')

subplot(4,4,[15,16])
plot(tHour(:,1), result(:,7), '.:'); xlim([1 length(tHour)]); grid on;
ylabel('Number of Satellites')

if length(GPGGA_sa) >= 1
    [length(GPGGA_dgps(:,1)), length(GPGGA_sa(:,1))]
else
    [length(GPGGA_dgps(:,1)), 0]
end