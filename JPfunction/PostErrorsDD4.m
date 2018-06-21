function [DDdXYZ, DDdNEV, DDdis, DDrms, DDstd] = PostErrorsDD4(estm, Base, A, B, ymin, ymax, SYS); 
%
% [DDdXYZ, DDdNEV, DDdis, DDrms, DDstd] = PostErrorsDD4(estm, Base, Truedis, ymin, ymax, visiSat); 
%
% <input>   estm: Rover XYZ [nx8]
%           base: Base XYZ [nx3]
%           Truedis: True distance between Base and Rover
%           ymin : limitation of subplot(3,2,1)
%           ymax : limitation of subplot(3,2,1)
%           SYS: Matrix of Satellites
%
% <output>  DDdXYZ: XYZ differences wrt Base
%           DDdNEV: NEV differences wrt Base
%           DDdis: 2D(1), 3D(2) distance between Base and Rover [n X 2]
%           DDrms: RMS of DDdNE_rms(1), DDdV_rms(2), DDd3_rms(3)
%           DDstd: Standard deviation of DDdNE_std(1), DDd3_std(2),
%
% Modified by JOON, 10/02/2017
%

%% 참값
Truedis = norm(B-A);
base = xyz2gd(A);
TrueNEV = xyz2topo((B-A), base(1), base(2));
% % TrueN = 
BsRvDis = Truedis;       % True Distance

FinalTTs = intersect(estm(:,1), Base(:,1));

for aa = 1:length(FinalTTs)
    base = Base(find(Base(:,1) == FinalTTs(aa)),2:4);
    Estm = estm(find(estm(:,1) == FinalTTs(aa)),2:4);
    Basegd = xyz2gd(base);
    DDdXYZ(aa,:) = Estm - base;
    DDdNEV(aa,:) = xyz2topo(DDdXYZ(aa,:), Basegd(1), Basegd(2));
    DDdN(aa,:) = DDdNEV(aa,1); 
    DDdE(aa,:) = DDdNEV(aa,2); 
    DDdV(aa,:) = DDdNEV(aa,3);
%     DDd3(aa,:) = (sqrt(DDdXYZ(aa,1)^2 + DDdXYZ(aa,2)^2 + DDdXYZ(aa,3)^2));
%     DDdNE(aa,:) = (sqrt(DDdXYZ(aa,1)^2 + DDdXYZ(aa,2)^2));
    DDdNE(aa,:) = (sqrt(DDdN(aa,1)^2 + DDdE(aa,1)^2));
    DDd3(aa,:) = (sqrt(DDdN(aa,1)^2 + DDdE(aa,1)^2 + DDdV(aa,1)^2));
    TrueDis(aa,1) = BsRvDis;
    
    
    
    
end
DDdis = [DDdNE, DDd3, TrueDis];
%% 각 성분별 RMS 계산
% [DDdNE_rms DDdV_rms DDd3_rms] = RMS(DDdNEV);
% [DDdNE_rms DDdV_rms DDd3_rms] = RMS(DDdXYZ);            
% [DDdNE_rms DDdV_rms DDd3_rms] = MEAN(DDdXYZ);
[DDdNE_rms DDdV_rms DDd3_rms] = MEAN(DDdNEV);
DDrms = [DDdNE_rms, DDdV_rms, DDd3_rms];
mean2D = zeros(length(estm),1); mean3D = zeros(length(estm),1);
DDdNE_std = std(DDdNE); DDd3_std = std(DDd3);
DDdN_std = std(DDdNEV(:,1)); DDdE_std = std(DDdNEV(:,2)); DDdV_std = std(DDdNEV(:,3));
DDstd = [DDdNE_std, DDd3_std, DDdN_std, DDdE_std, DDdV_std];

tHour = mod(FinalTTs, 86400)/3600;

%% Plot

figure(199)
subplot(3,2,1)
hold on; grid on;
% xlim([0 length(estm)])
xlim([min(tHour) max(tHour)])
ylim([ymin ymax])
plot(tHour, DDdNE,'b.:');
% plot(DDd3,'b.:');
plot(TrueDis,'r-');
drawnow
xlabel({['Double Differencing  2D = ', num2str(decimal(DDdNE_rms)), '   std =', num2str(decimal(DDdNE_std))]}); 
% xlabel({['Double Differencing  2D = ', num2str(DDdNE_rms), '   std =', num2str(DDdNE_std)]}); 
% xlabel({['Double Differencing  3D = ', num2str(DDd3_rms), '   std =', num2str(DDd3_std)]}); 
ylabel('2D Distance(meter)');
legend('2D-dis','TrueDis','Location','Best')
% legend('3D-Dis','TrueDis')

subplot(3,2,3)
hold on; grid on;
% xlim([0 length(estm)])
xlim([min(tHour) max(tHour)])
ylim([ymin ymax])
% plot(DDdNE,'b.:');
plot(tHour, DDd3,'b.:');
plot(TrueDis,'r-');
drawnow
xlabel({['Double Differencing  3D = ', num2str(decimal(DDd3_rms)), '   std =', num2str(decimal(DDd3_std))]}); 
% xlabel({['Double Differencing  2D = ', num2str(decimal(DDdNE_rms)), '   std =', num2str(decimal(DDdNE_std))],...
%         ['Double Differencing  3D = ', num2str(decimal(DDd3_rms)), '   std =', num2str(decimal(DDd3_std))]}); 
% xlabel({['Double Differencing  2D = ', num2str(DDdNE_rms), '   std =', num2str(DDdNE_std)]}); 
% xlabel({['Double Differencing  3D = ', num2str(DDd3_rms), '   std =', num2str(DDd3_std)]}); 
ylabel('3D Distance(meter)');
legend('3D-dis','TrueDis','Location','Best')
% legend('3D-Dis','TrueDis')

subplot(3,2,2)
hold on; grid on;
% xlim([0 length(estm)])
xlim([min(tHour) max(tHour)])
% ylim([-5 5])
% plot(DDdXYZ(:,1),'r.:');
plot(tHour, DDdNEV(:,1),'b.:');
plot([min(tHour) max(tHour)],[TrueNEV(1),TrueNEV(1)],'r-')
xlabel({['Epoch'], ['Standard Deviation(dN) = ',num2str(decimal(DDdN_std))]});
% ylabel('dX(meter)','fontsize',15);
ylabel('dN(meter)','fontsize',15);
% legend('dX')
legend('dN','Location','Best')

subplot(3,2,4)
hold on; grid on;
% xlim([0 length(estm)])
xlim([min(tHour) max(tHour)])
% ylim([-5 5])
% plot(DDdXYZ(:,2),'b.:');
plot(tHour, DDdNEV(:,2),'b.:');
plot([min(tHour) max(tHour)],[TrueNEV(2),TrueNEV(2)],'r-')
xlabel({['Epoch'], ['Standard Deviation(dE) = ',num2str(decimal(DDdE_std))]});
% ylabel('dY(meter)','fontsize',15);
ylabel('dE(meter)','fontsize',15);
% legend('dY')
legend('dE','Location','Best')

subplot(3,2,5)
hold on; grid on;
% xlim([0 length(estm)])
xlim([min(tHour) max(tHour)])
if SYS == 1
    stairs(tHour, estm(:,7),'b.:');
    legend('GPS','Location','Best')
elseif SYS == 2
    stairs(tHour, estm(:,8),'r.:');
    legend('BDS','Location','Best')
else
    stairs(tHour, estm(:,7),'b.:');
    stairs(tHour, estm(:,8),'r.:');
    legend('GPS','BDS','Location','Best')
end
% stairs(tHour, visiSat(:,8),'b.:');
% plot(DDdXYZ(:,1),'r.:');
% plot(DDdXYZ(:,2),'b.:');
% plot(DDdXYZ(:,3),'k.:');
% ylim([-5 5]);
% ylim([8 20]);
xlabel('Epoch','fontsize',15);
% ylabel('meter','fontsize',15);
ylabel('nums of sats','fontsize',15);
% legend('dX','dY','dZ')


subplot(3,2,6)
hold on; grid on;
% xlim([0 length(estm)])
xlim([min(tHour) max(tHour)])
ylim([-5 5])
% plot(DDdXYZ(:,3),'k.:');
plot(tHour, DDdNEV(:,3),'b.:');
plot([min(tHour) max(tHour)],[TrueNEV(3),TrueNEV(3)],'r-')
xlabel({['Epoch'], ['Standard Deviation(dV) = ',num2str(decimal(DDdV_std))]});
% ylabel('dZ(meter)','fontsize',15);
ylabel('dV(meter)','fontsize',15);
% legend('dZ')
legend('dV')