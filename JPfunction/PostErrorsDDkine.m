function [DDdXYZ, DDdXYZ_vrs, DDdNEV, DDdNEV_vrs, DDdis, DDrms, result] = PostErrorsDDkine(estm, Base, Base_vrs, Rover_vrs, ymin, ymax, SYS); 
%
% [DDdXYZ, DDdNEV, DDdis, DDrms, DDstd] = PostErrorsDDkine(estm, Base, Truedis, ymin, ymax, SYS); 
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
% clear all
% load('posterrorkinetest.mat')
% ymin = -10; ymax = 10;
%% 참값
Base_vrs(:,1) = Base_vrs(:,1) + 18;
Rover_vrs(:,1) = Rover_vrs(:,1) + 18 ;
FinalTTs1 = intersect(Base_vrs(:,1), Rover_vrs(:,1));
FinalTTs2 = intersect(estm(:,1), Base(:,1));
FinalTTs = intersect(FinalTTs1, FinalTTs2);
for aa = 1:length(FinalTTs)
    gs = FinalTTs(aa);
    base = Base(find(Base(:,1) == FinalTTs(aa)),2:4);
    base_vrs = Base_vrs(find(Base_vrs(:,1) == FinalTTs(aa)),5:7);
    Estm = estm(find(estm(:,1) == FinalTTs(aa)),2:4);
    estout(aa,:) = estm(find(estm(:,1) == FinalTTs(aa)),:);
    rover_vrs = Rover_vrs(find(Rover_vrs(:,1) == FinalTTs(aa)),5:7);
    Basegd = xyz2gd(base);
    Base_vrs_gd = Base_vrs(find(Base_vrs(:,1) == FinalTTs(aa)),2:4);
    GD(aa,1:4) = [gs Basegd];
    DDdXYZ(aa,:) = Estm - base;
    DDdXYZ_vrs(aa,:) = rover_vrs(1,:) - base_vrs(1,:);
    DDdNEV(aa,:) = xyz2topo(DDdXYZ(aa,:), Basegd(1), Basegd(2));
    DDdNEV_vrs(aa,:) = xyz2topo(DDdXYZ_vrs(aa,:), Base_vrs_gd(1), Base_vrs_gd(2));
    NEV(aa,1:4) = [gs DDdNEV(aa,:)];
    DDdN(aa,:) = DDdNEV(aa,1); 
    DDdN_vrs(aa,:) = DDdNEV_vrs(aa,1); 
    DDdE(aa,:) = DDdNEV(aa,2); 
    DDdE_vrs(aa,:) = DDdNEV_vrs(aa,2); 
    DDdV(aa,:) = DDdNEV(aa,3);
    DDdV_vrs(aa,:) = DDdNEV_vrs(aa,3);
    DDdNE(aa,:) = [gs, (sqrt(DDdN(aa,1)^2 + DDdE(aa,1)^2))];
    DDdNE_vrs(aa,:) = [gs, (sqrt(DDdN_vrs(aa,1)^2 + DDdE_vrs(aa,1)^2))];
    DDd3(aa,:) = (sqrt(DDdN(aa,1)^2 + DDdE(aa,1)^2 + DDdV(aa,1)^2));
    DDd3_vrs(aa,:) = (sqrt(DDdN_vrs(aa,1)^2 + DDdE_vrs(aa,1)^2 + DDdV_vrs(aa,1)^2));
%     DDd3(aa,:) = (sqrt(DDdXYZ(aa,1)^2 + DDdXYZ(aa,1)^2 + DDdXYZ(aa,1)^2));
%     DDd3_vrs(aa,:) = (sqrt(DDdXYZ_vrs(aa,1)^2 + DDdXYZ_vrs(aa,1)^2 + DDdXYZ_vrs(aa,1)^2));
    TrueDis(aa,1:3) = [gs, DDdNE_vrs(aa,2), DDd3_vrs(aa,1)];
    result(aa,1:3) = [gs, DDdNE_vrs(aa,2) - DDdNE(aa,2), DDd3_vrs(aa,1) - DDd3(aa,1)];
    
end
estm =estout;
DDdis = [TrueDis, DDdNE(:,2), DDd3];
DDdis_2D(:,1) = smooth(DDdis(:,2)-DDdis(:,4),0.1,'loess');
%% 각 성분별 RMS 계산
% [DDdNE_rms DDdV_rms DDd3_rms] = RMS(DDdNEV-DDdNEV_vrs);
% [DDdNE_rms DDdV_rms DDd3_rms] = RMS(DDdXYZ);            
% [DDdNE_rms DDdV_rms DDd3_rms] = MEAN(DDdXYZ);
[DDdNE_rms DDdV_rms DDd3_rms] = MEAN((DDdNEV-DDdNEV_vrs));
DDrms = [rms(result(:,2)), rms(result(:,3))];
% DDdNE_std = std(result(:,5)); DDd3_std = std(result(:,6));
% DDdN_std = std(result(:,2)); DDdE_std = std(result(:,3)); DDdV_std = std(result(:,4));
% DDstd = [DDdNE_std, DDd3_std, DDdN_std, DDdE_std, DDdV_std];

tHour = mod(FinalTTs, 86400)/3600;

%% Plot

figure()
subplot(3,1,1)
hold on; grid on;
% xlim([0 length(estm)])
xlim([min(tHour) max(tHour)])
ylim([ymin ymax])
plot(tHour, DDdis(:,2)-DDdis(:,4),'b.:');
% plot(tHour, DDdis_2D(:,1),'b.:');
drawnow
xlabel({['Double Differencing  2D = ', num2str(decimal(rms(DDdis(:,2)-DDdis(:,4)))), '   std =', num2str(decimal(std(DDdis(:,2)-DDdis(:,4))))]}); 

ylabel('2D Distance(meter)');
legend('2D-dis')
% legend('3D-Dis','TrueDis')

subplot(3,1,2)
hold on; grid on;
% xlim([0 length(estm)])
xlim([min(tHour) max(tHour)])
ylim([ymin ymax])
% plot(DDdNE,'b.:');
plot(tHour, DDdis(:,3)-DDdis(:,5),'b.:');
% plot(TrueDis,'k-');
drawnow
xlabel({['Double Differencing  3D = ', num2str(decimal(rms(DDdis(:,3)-DDdis(:,5)))), '   std =', num2str(decimal(std(DDdis(:,3)-DDdis(:,5))))]}); 
ylabel('3D Distance(meter)');
legend('3D-dis')
% legend('3D-Dis','TrueDis')

% subplot(3,2,2)
% hold on; grid on;
% % xlim([0 length(estm)])
% xlim([min(tHour) max(tHour)])
% % ylim([-5 5])
% % plot(DDdXYZ(:,1),'r.:');
% plot(tHour, result(:,2),'b.:');
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
% plot(tHour, result(:,3),'b.:');
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
    stairs(tHour, estm(:,8),'r.:');
    stairs(tHour, estm(:,9),'g.:');
    legend('GPS','BDS')
elseif SYS == 6
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
% plot(tHour, result(:,4),'b.:');
% xlabel({['Epoch'], ['Standard Deviation(dV) = ',num2str(decimal(DDdV_std))]});
% % ylabel('dZ(meter)','fontsize',15);
% ylabel('dV(meter)','fontsize',15);
% % legend('dZ')
% legend('dV')