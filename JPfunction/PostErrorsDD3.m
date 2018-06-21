function [DDdXYZ, DDdNEV, DDdis, DDrms, DDstd] = PostErrorsDD3(estm, Base, TruePosBs, TruePosRv); 
%
% function [DDdXYZ, DDdNEV, DDdis, rms, std] = PostErrorsDD3(estm, Base, TruePosBs, TruePosRv);  
%
% <input>   estm: Rover XYZ [nx3]
%           base: Base XYZ [nx3]
%           TruePosBs: Base True Position XYZ [1x3]
%           TruePosRv: Rover True Position XYZ [1x3]
%
% <output>  DDdXYZ: XYZ differences wrt Base
%           DDdNEV: NEV differences wrt Base
%           DDdis: 2D(1), 3D(2) distance between Base and Rover [n X 2]
%           DDrms: RMS of DDdNE_rms(1), DDdV_rms(2), DDd3_rms(3)
%           DDstd: Standard deviation of DDdNE_std(1), DDd3_std(2),
%           DDdN_std(3), DDdE_std(4), DDdV_std(5)
%
% Modified by JOON, 24/02/2016
%

%% 참값
BsRvDis = sqrt((TruePosBs(1)-TruePosRv(1))^2+...
    (TruePosBs(2)-TruePosRv(2))^2+(TruePosBs(3)-TruePosRv(3))^2);       % True Distance
ddxyz = TruePosBs -TruePosRv

for aa = 1:length(estm)
    Basegd = xyz2gd(Base(aa,2:4));
    DDdXYZ(aa,:) = estm(aa,2:4) - TruePosRv;%Base(aa,2:4);
    DDdNEV(aa,:) = xyz2topo(DDdXYZ(aa,:), Basegd(1), Basegd(2));
    DDdN(aa,:) = DDdNEV(aa,1); 
    DDdE(aa,:) = DDdNEV(aa,2); 
    DDdV(aa,:) = DDdNEV(aa,3);
    DDdNE(aa,:) = (sqrt(DDdXYZ(aa,1)^2 + DDdXYZ(aa,2)^2));
    DDd3(aa,:) = (sqrt(DDdXYZ(aa,1)^2 + DDdXYZ(aa,2)^2 + DDdXYZ(aa,3)^2));
    TrueDis(aa,1) = BsRvDis;
    TruedX(aa,1) = ddxyz(1);
    TruedY(aa,1) = ddxyz(2);
    TruedZ(aa,1) = ddxyz(3);
     
    
end
DDdis = [DDdNE, DDd3, TrueDis];
%% 각 성분별 RMS 계산
[DDdNE_rms DDdV_rms DDd3_rms] = RMS(DDdNEV);
DDrms = [DDdNE_rms, DDdV_rms, DDd3_rms];
mean2D = zeros(length(estm),1); mean3D = zeros(length(estm),1);
DDdNE_std = std(DDdNE); DDd3_std = std(DDd3);
DDdN_std = std(DDdNEV(:,1)); DDdE_std = std(DDdNEV(:,2)); DDdV_std = std(DDdNEV(:,3));
DDstd = [DDdNE_std, DDd3_std, DDdN_std, DDdE_std, DDdV_std];

%% Plot

figure(99)
subplot(3,2,3)
hold on; grid on;
xlim([0 length(estm)])
ylim([0 5])
% plot(DDdNE,'r.:', 'MarkerSize',10);
plot(DDd3,'b.:', 'MarkerSize',10);
plot(TrueDis,'r-');
drawnow
xlabel('Hours');
ylabel('Baseline Length');
% legend('2D(dNE) Distance','3D(d3D) Distance','True Distance')

% subplot(3,2,1)
% hold on; grid on;
% xlim([0 length(estm)])
% plot(DDdN,'r.:');
% plot(DDdE,'b.:');
% plot(DDdV,'k.:');
% xlabel('Epoch','fontsize',15);
% ylabel('meter','fontsize',15);
% legend('dN','dE','dV')

subplot(3,2,2)
hold on; grid on;
xlim([0 length(estm)])
ylim([-5 5])
plot(DDdXYZ(:,1),'b.:', 'MarkerSize',10);
plot(TruedX,'r-');
% xlabel({['Epoch'], ['Standard Deviation(dN) = ',num2str(DDdN_std)]});
ylabel('dX(meter)','fontsize',15);
% legend('dX')

subplot(3,2,4)
hold on; grid on;
xlim([0 length(estm)])
ylim([-5 5])
plot(DDdXYZ(:,2),'b.:', 'MarkerSize',10);
plot(TruedY,'r-');
% xlabel({['Epoch'], ['Standard Deviation(dE) = ',num2str(DDdE_std)]});
ylabel('dY(meter)','fontsize',15);
% legend('dY')

subplot(3,2,6)
hold on; grid on;
xlim([0 length(estm)])
ylim([-5 5])
plot(DDdXYZ(:,3),'b.:', 'MarkerSize',10);
plot(TruedZ,'r-');
% xlabel({['Epoch'], ['Standard Deviation(dV) = ',num2str(DDdV_std)]});
ylabel('dZ(meter)','fontsize',15);
% legend('dZ')