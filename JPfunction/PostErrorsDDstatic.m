function [DDdXYZ, DDdNEV, DDdis, DDrms, DDstd] = PostErrorsDDstatic(estm, Base, Truedis, ymin, ymax, SYS); 
%
% [DDdXYZ, DDdNEV, DDdis, DDrms, DDstd] = PostErrorsDDstatic(estm, Base, Truedis, ymin, ymax, SYS); 
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
BsRvDis = Truedis;       % True Distance
estm(:,1) = round(estm(:,1));
Base(:,1) = round(Base(:,1));
FinalTTs = intersect(estm(:,1), Base(:,1));

for aa = 1:length(FinalTTs)
    gs = FinalTTs(aa);
    base = Base(find(Base(:,1) == FinalTTs(aa)),2:4);
    Estm = estm(find(estm(:,1) == FinalTTs(aa)),2:4);
    estout(aa,:) = estm(find(estm(:,1) == FinalTTs(aa)),:);
    Basegd = xyz2gd(base);
    GD(aa,1:4) = [gs Basegd];
    DDdXYZ(aa,:) = Estm - base;
    DDdNEV(aa,:) = xyz2topo(DDdXYZ(aa,:), Basegd(1), Basegd(2));
    NEV(aa,1:4) = [gs DDdNEV(aa,:)];
    DDdN(aa,:) = DDdNEV(aa,1); 
    DDdE(aa,:) = DDdNEV(aa,2); 
    DDdV(aa,:) = DDdNEV(aa,3);
%     DDd3(aa,:) = (sqrt(DDdXYZ(aa,1)^2 + DDdXYZ(aa,2)^2 + DDdXYZ(aa,3)^2));
%     DDdNE(aa,:) = (sqrt(DDdXYZ(aa,1)^2 + DDdXYZ(aa,2)^2));
    DDdNE(aa,:) = [gs, (sqrt(DDdN(aa,1)^2 + DDdE(aa,1)^2))];
    DDd3(aa,:) = (sqrt(DDdN(aa,1)^2 + DDdE(aa,1)^2 + DDdV(aa,1)^2));
    TrueDis(aa,1) = BsRvDis;
    
    
    
    
end
estm = estout;
DDdis = [DDdNE, DDd3, TrueDis];
%% 각 성분별 RMS 계산
[DDdNE_rms DDdV_rms DDd3_rms] = RMS(DDdNEV);
% [DDdNE_rms DDdV_rms DDd3_rms] = RMS(DDdXYZ);            
% [DDdNE_rms DDdV_rms DDd3_rms] = MEAN(DDdXYZ);
% [DDdNE_rms DDdV_rms DDd3_rms] = MEAN(DDdNEV);
DDrms = [DDdNE_rms, DDdV_rms, DDd3_rms];
mean2D = zeros(length(estm),1); mean3D = zeros(length(estm),1);
DDdNE_std = std(DDdNE); DDd3_std = std(DDd3);
DDdN_std = std(DDdNEV(:,1)); DDdE_std = std(DDdNEV(:,2)); DDdV_std = std(DDdNEV(:,3));
DDstd = [DDdNE_std, DDd3_std, DDdN_std, DDdE_std, DDdV_std];

tHour = mod(FinalTTs, 86400)/3600;

%% Plot

figure()
subplot(3,1,1)
hold on; grid on;
% xlim([0 length(estm)])
xlim([min(tHour) max(tHour)])
ylim([ymin ymax])
plot(tHour, DDdNE,'b.:');
% plot(DDd3,'b.:');
plot(TrueDis,'k-');
drawnow
xlabel({['Double Differencing  2D = ', num2str(decimal(DDdNE_rms)), '   std =', num2str(decimal(DDdNE_std(:,2)))]}); 
% xlabel({['Double Differencing  2D = ', num2str(DDdNE_rms), '   std =', num2str(DDdNE_std)]}); 
% xlabel({['Double Differencing  3D = ', num2str(DDd3_rms), '   std =', num2str(DDd3_std)]}); 
ylabel('2D Distance(meter)');
legend('2D-dis','TrueDis')
% legend('3D-Dis','TrueDis')

subplot(3,1,2)
hold on; grid on;
% xlim([0 length(estm)])
xlim([min(tHour) max(tHour)])
ylim([ymin ymax])
% plot(DDdNE,'b.:');
plot(tHour, DDd3,'b.:');
plot(TrueDis,'k-');
drawnow
xlabel({['Double Differencing  3D = ', num2str(decimal(DDd3_rms)), '   std =', num2str(decimal(DDd3_std))]}); 
% xlabel({['Double Differencing  2D = ', num2str(decimal(DDdNE_rms)), '   std =', num2str(decimal(DDdNE_std))],...
%         ['Double Differencing  3D = ', num2str(decimal(DDd3_rms)), '   std =', num2str(decimal(DDd3_std))]}); 
% xlabel({['Double Differencing  2D = ', num2str(DDdNE_rms), '   std =', num2str(DDdNE_std)]}); 
% xlabel({['Double Differencing  3D = ', num2str(DDd3_rms), '   std =', num2str(DDd3_std)]}); 
ylabel('3D Distance(meter)');
legend('3D-dis','TrueDis')
% legend('3D-Dis','TrueDis')

% subplot(3,2,2)
% hold on; grid on;
% % xlim([0 length(estm)])
% xlim([min(tHour) max(tHour)])
% % ylim([-5 5])
% % plot(DDdXYZ(:,1),'r.:');
% plot(tHour, DDdNEV(:,1),'b.:');
% xlabel({['Epoch'], ['Standard Deviation(dN) = ',num2str(decimal(DDdN_std))]});
% % ylabel('dX(meter)','fontsize',15);
% ylabel('dN(meter)','fontsize',15);
% % legend('dX')
% legend('dN')

% subplot(3,2,4)
% hold on; grid on;
% % xlim([0 length(estm)])
% xlim([min(tHour) max(tHour)])
% % ylim([-5 5])
% % plot(DDdXYZ(:,2),'b.:');
% plot(tHour, DDdNEV(:,2),'b.:');
% xlabel({['Epoch'], ['Standard Deviation(dE) = ',num2str(decimal(DDdE_std))]});
% % ylabel('dY(meter)','fontsize',15);
% ylabel('dE(meter)','fontsize',15);
% % legend('dY')
% legend('dE')

subplot(3,1,3)
hold on; grid on;
% xlim([0 length(estm)])
xlim([min(tHour) max(tHour)])
if SYS == 1
    stairs(tHour, estm(:,7),'b.:');
    legend('GPS')
elseif SYS == 2
    stairs(tHour, estm(:,8),'r.:');
    legend('BDS')
elseif SYS == 3
    stairs(tHour, estm(:,9),'g.:');
    legend('GLO')
elseif SYS == 4
    stairs(tHour, estm(:,7),'b.:');
    stairs(tHour, estm(:,8),'r.:');
    legend('GPS','BDS')
elseif SYS == 5
    stairs(tHour, estm(:,7),'b.:');
    stairs(tHour, estm(:,9),'g.:');
    legend('GPS','GLO')
elseif SYS == 6
    stairs(tHour, estm(:,8),'r.:');
    stairs(tHour, estm(:,9),'g.:');
    legend('BDS','GLO')
elseif SYS == 7
    stairs(tHour, estm(:,7),'b.:');
    stairs(tHour, estm(:,8),'r.:');
    stairs(tHour, estm(:,9),'g.:');
    legend('GPS','BDS', 'GLO')
end

xlabel('Epoch','fontsize',15);
ylabel('nums of sats','fontsize',15);



% subplot(3,2,6)
% hold on; grid on;
% % xlim([0 length(estm)])
% xlim([min(tHour) max(tHour)])
% ylim([-5 5])
% % plot(DDdXYZ(:,3),'k.:');
% plot(tHour, DDdNEV(:,3),'b.:');
% xlabel({['Epoch'], ['Standard Deviation(dV) = ',num2str(decimal(DDdV_std))]});
% % ylabel('dZ(meter)','fontsize',15);
% ylabel('dV(meter)','fontsize',15);
% % legend('dZ')
% legend('dV')