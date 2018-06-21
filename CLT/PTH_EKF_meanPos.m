clear all; close all;
tic;

%% 좌우 안테나 위치 결과 load
file_L = 'test3_L_EstOut.mat';
file_R = 'test3_R_EstOut.mat';
load(file_L); load(file_R);

if file_L(end-5:end-4) == 'el'
    EstOut_L = EstOut_L_el;
    EstOut_R = EstOut_R_el;
    title_txt = 'EKF with el weighting';
else
    title_txt = 'EKF with SNR weighting';
end
%% VRS load
load('PTCO1_170116_adm.txt');
vrs= PTCO1_170116_adm;
vrs(:,1) = vrs(:,1)+18;

%% 기준 좌표
TruePos = [-3058799.61420451,4083265.35912516,3814946.87192938];
gd = xyz2gd(TruePos); TrueLat = gd(1); TrueLon = gd(2);

%% 공통시간
FinalTTs = intersect(EstOut_L(:, 1), EstOut_R(:, 1));

for i=1:length(FinalTTs)
    gs = FinalTTs(i);
    vrs_coordi = vrs(find(vrs(:,1)==gs),5:7);
    left = EstOut_L(find(EstOut_L(:,1) == gs),2:4);
    right = EstOut_R(find(EstOut_R(:,1) == gs),2:4);
    pos = (left + right)/2;
    
    %% dNEV 계산
    dXYZ_vrs = vrs_coordi - TruePos;
    dXYZ_pos = pos - TruePos;
    dXYZ_pos_vrs = pos - vrs_coordi;
    dNEV_vrs = xyz2topo(dXYZ_vrs, TrueLat, TrueLon);
    dNEV_pos = xyz2topo(dXYZ_pos, TrueLat, TrueLon);
    dNEV_pos_vrs = xyz2topo(dXYZ_pos_vrs, TrueLat, TrueLon);
    dN_vrs = dNEV_vrs(:,1); dE_vrs = dNEV_vrs(:,2); dV_vrs = dNEV_vrs(:,3);
    dN_pos = dNEV_pos(:,1); dE_pos = dNEV_pos(:,2); dV_pos = dNEV_pos(:,3);
    dN_pos_vrs = dNEV_pos_vrs(:,1); dE_pos_vrs = dNEV_pos_vrs(:,2); dV_pos_vrs = dNEV_pos_vrs(:,3);
    dNE_vrs = sqrt(dN_vrs^2 + dE_vrs^2);
    dNE_pos = sqrt(dN_pos^2 + dE_pos^2);
    dNE_pos_vrs = sqrt(dN_pos_vrs^2 + dE_pos_vrs^2);
    d2D_vrs(i,1) = dNE_vrs;
    d2D_pos(i,1) = dNE_pos;
    d2D_pos_vrs(i,1) = dNE_pos_vrs;
    d3_vrs = sqrt(dN_vrs.^2 + dE_vrs.^2 + dV_vrs.^2);
    d3_pos = sqrt(dN_pos.^2 + dE_pos.^2 + dV_pos.^2);
    d3_pos_vrs = sqrt(dN_pos_vrs.^2 + dE_pos_vrs.^2 + dV_pos_vrs.^2);
    d3D_vrs(i,1) = d3_vrs;
    d3D_pos(i,1) = d3_pos;
    d3D_pos_vrs(i,1) = d3_pos_vrs;
    result_vrs(i,1:6) = [gs, dN_vrs, dE_vrs, dV_vrs, dNE_vrs, d3_vrs];
    result_pos(i,1:11) = [gs, dN_pos, dE_pos, dV_pos, dNE_pos, d3_pos, dN_pos_vrs, dE_pos_vrs, dV_pos_vrs, dNE_pos_vrs, d3_pos_vrs];
end

%% 실험 세트 분리
vrs1 = result_vrs(1:592,:);
vrs2 = result_vrs(620:1140,:);
vrs3 = result_vrs(1158:1636,:);
pos1 = result_pos(1:592,:);
pos2 = result_pos(620:1140,:);
pos3 = result_pos(1158:1636,:);

figure(1)
suptitle(title_txt)
tHour = mod(result_pos(:,1), 86400);
tHour = tHour/3600;

% subplot(3,5,[1,2,3,6,7,8,11,12,13])
plot(result_pos(:,8),result_pos(:,7),'bo')
grid on; hold on;
plot([-20,20],[0,0],'r-')
plot([0,0],[-20,20],'r-')
axis([-5 5 -5 5])
axis square
xlabel(['H RMSE = ', num2str(decimal(rms(result_pos(:,10))))]);

% subplot(3,5,[4,5])
% plot(tHour,result_pos(:,7),'b.:')
% grid on; hold on;
% xlim([tHour(1) tHour(length(tHour))]);
% ylim([-5 5])
% ylabel('\Delta N')
% xlabel('Hour')
% 
% subplot(3,5,[9,10])
% plot(tHour,result_pos(:,8),'b.:')
% grid on; hold on;
% xlim([tHour(1) tHour(length(tHour))]);
% ylim([-5 5])
% ylabel('\Delta E')
% xlabel('Hour')
% 
% subplot(3,5,[14,15])
% plot(tHour,result_pos(:,9),'b.:')
% grid on; hold on;
% xlim([tHour(1) tHour(length(tHour))]);
% ylim([-5 5])
% ylabel('\Delta V')
% xlabel('Hour')
% tHour = 0;

figure(2)
suptitle(title_txt)
plot(result_pos(:,3),result_pos(:,2),'bo')
hold on; grid on;
plot(result_vrs(:,3),result_vrs(:,2),'ro')
axis([-60 -25 30 58])
axis equal
legend('PPSoln','VRS')

%% 정지시점 결과 추출
cnt=1;
for i = 2:length(vrs1(:,1))
    gs = vrs1(i-1,1);
    temp2 = vrs1(i,5);
    temp1 = vrs1(i-1,5);
    diff = abs(temp1 -temp2);
    if diff < 0.02
        vrs1_stop(cnt,:) = vrs1(i-1,:);
        pos1_stop(cnt,:) = pos1(find(pos1(:,1) == gs),:);
        cnt = cnt + 1;
    end
end
cnt=1;
for i = 2:length(vrs2(:,1))
    gs = vrs2(i-1,1);
    temp2 = vrs2(i,5);
    temp1 = vrs2(i-1,5);
    diff = abs(temp1 -temp2);
    if diff < 0.02
        vrs2_stop(cnt,:) = vrs2(i-1,:);
        pos2_stop(cnt,:) = pos2(find(pos2(:,1) == gs),:);
        cnt = cnt + 1;
    end
end
cnt=1;
for i = 2:length(vrs3(:,1))
    gs = vrs3(i-1,1);
    temp2 = vrs3(i,5);
    temp1 = vrs3(i-1,5);
    diff = abs(temp1 -temp2);
    if diff < 0.02
        vrs3_stop(cnt,:) = vrs3(i-1,:);
        pos3_stop(cnt,:) = pos3(find(pos3(:,1) == gs),:);
        cnt = cnt + 1;
    end
end

VRS_STOP = [vrs1_stop;vrs2_stop;vrs3_stop];
POS_STOP = [pos1_stop;pos2_stop;pos3_stop];


figure(3)
suptitle([title_txt,' STOP'])
tHour = mod(POS_STOP(:,1), 86400);
tHour = tHour/3600;

% subplot(3,5,[1,2,3,6,7,8,11,12,13])
plot(POS_STOP(:,8),POS_STOP(:,7),'bo')
grid on; hold on;
plot([-20,20],[0,0],'r-')
plot([0,0],[-20,20],'r-')
axis([-5 5 -5 5])
axis square
xlabel(['H RMSE = ', num2str(decimal(rms(POS_STOP(:,10))))]);

figure(4)
suptitle([title_txt,' STOP'])
hold on; grid on;
plot(vrs2(:,3), vrs2(:,2),'r:','LineWidth',0.5);
plot(VRS_STOP(:,3),VRS_STOP(:,2),'ro')
plot(POS_STOP(:,3),POS_STOP(:,2),'bo')
axis([-60 -25 30 58])
axis equal
legend('PPSoln','VRS')

figure(5)
subplot(3,1,1)
title([title_txt, '(set 1)'])
hold on; grid on;
plot(vrs1(:,3),vrs1(:,2),'r:','LineWidth',0.5)
plot(vrs1_stop(:,3), vrs1_stop(:,2),'r.','MarkerSize',10);
plot(pos1_stop(:,3), pos1_stop(:,2),'b.','MarkerSize',10);
legend('VRS Path','VRS(Stop)','Pos(Stop)','Location','Best')
xlabel({['dE'],['2D RMSE = ',num2str(decimal(rms(pos1_stop(:,10)))),' m']})
ylabel('dN')
axis equal

subplot(3,1,2)
title([title_txt, '(set 2)'])
hold on; grid on;
plot(vrs2(:,3),vrs2(:,2),'r:','LineWidth',0.5)
plot(vrs2_stop(:,3), vrs2_stop(:,2),'r.','MarkerSize',10);
plot(pos2_stop(:,3), pos2_stop(:,2),'b.','MarkerSize',10);
legend('VRS Path','VRS(Stop)','Pos(Stop)','Location','Best')
xlabel({['dE'],['2D RMSE = ',num2str(decimal(rms(pos2_stop(:,10)))),' m']})
ylabel('dN')
axis equal

subplot(3,1,3)
title([title_txt, '(set 3)'])
hold on; grid on;
plot(vrs2(:,3),vrs2(:,2),'r:','LineWidth',0.5)
plot(vrs3_stop(:,3), vrs3_stop(:,2),'r.','MarkerSize',10);
plot(pos3_stop(:,3), pos3_stop(:,2),'b.','MarkerSize',10);
legend('VRS Path','VRS(Stop)','Pos(Stop)','Location','Best')
xlabel({['dE'],['2D RMSE = ',num2str(decimal(rms(pos3_stop(:,10)))),' m']})
ylabel('dN')
axis equal

sigma1 = decimal(sigma_range(result_pos(:,10),1))
sigma2 = decimal(sigma_range(result_pos(:,10),2))
under1m = decimal((length(result_pos(find(result_pos(:,10) <= 1),10))/length(result_pos(:,10))) * 100)
under1_2m = decimal((length(result_pos(find(result_pos(:,10) <= 1.2),10))/length(result_pos(:,10))) * 100)

sigma1 = decimal(sigma_range(POS_STOP(:,10),1))
sigma2 = decimal(sigma_range(POS_STOP(:,10),2))
under1m = decimal((length(POS_STOP(find(POS_STOP(:,10) <= 1),10))/length(POS_STOP(:,10))) * 100)
under1_2m = decimal((length(POS_STOP(find(POS_STOP(:,10) <= 1.2),10))/length(POS_STOP(:,10))) * 100)


figure(6)
title([title_txt,' (set1)'])
hold on; grid on;
plot(vrs1(:,3), vrs1(:,2),'r.','MarkerSize',10);
plot(pos1(:,3), pos1(:,2),'b.','MarkerSize',10);
xlabel({['dE'],['2D RMSE = ',num2str(decimal(rms(pos1(:,10)))),' m']})
ylabel('dN')
axis equal
figure(7)
title([title_txt,' (set2)'])
hold on; grid on;
plot(vrs2(:,3), vrs2(:,2),'r.','MarkerSize',10);
plot(pos2(:,3), pos2(:,2),'b.','MarkerSize',10);
xlabel({['dE'],['2D RMSE = ',num2str(decimal(rms(pos2(:,10)))),' m']})
ylabel('dN')
axis equal
figure(8)
title([title_txt,' (set3)'])
hold on; grid on;
plot(vrs3(:,3), vrs3(:,2),'r.','MarkerSize',10);
plot(pos3(:,3), pos3(:,2),'b.','MarkerSize',10);
xlabel({['dE'],['2D RMSE = ',num2str(decimal(rms(pos3(:,10)))),' m']})
ylabel('dN')
axis equal

figure(9)
title([title_txt,' STOP(set1)'])
hold on; grid on;
plot(vrs1(:,3), vrs1(:,2),'r:','LineWidth',0.5);
plot(vrs1_stop(:,3), vrs1_stop(:,2),'r.','MarkerSize',10);
plot(pos1_stop(:,3), pos1_stop(:,2),'b.','MarkerSize',10);
xlabel({['dE'],['2D RMSE = ',num2str(decimal(rms(pos1_stop(:,10)))),' m']})
ylabel('dN')
axis equal
figure(10)
title([title_txt,' STOP(set2)'])
hold on; grid on;
plot(vrs2(:,3), vrs2(:,2),'r:','LineWidth',0.5);
plot(vrs2_stop(:,3), vrs2_stop(:,2),'r.','MarkerSize',10);
plot(pos2_stop(:,3), pos2_stop(:,2),'b.','MarkerSize',10);
xlabel({['dE'],['2D RMSE = ',num2str(decimal(rms(pos2_stop(:,10)))),' m']})
ylabel('dN')
axis equal
figure(11)
title([title_txt,' STOP(set3)'])
hold on; grid on;
plot(vrs3(:,3), vrs3(:,2),'r:','LineWidth',0.5);
plot(vrs3_stop(:,3), vrs3_stop(:,2),'r.','MarkerSize',10);
plot(pos3_stop(:,3), pos3_stop(:,2),'b.','MarkerSize',10);
xlabel({['dE'],['2D RMSE = ',num2str(decimal(rms(pos3_stop(:,10)))),' m']})
ylabel('dN')
axis equal