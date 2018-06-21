function [DDdXYZ, DDdNEV, DDdis, DDrms, DDstd] = PostErrorsDD(estm, Base, TruePosBs, TruePosRv); 
%
% function [DDdXYZ, DDdNEV, DDdis, rms, std] = PostErrorsDD(estm, Base, TruePosBs, TruePosRv);  
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
Truegd = xyz2gd(TruePosRv);
dBsRv = TruePosRv - TruePosBs;
TruedNEV = xyz2topo(dBsRv, Truegd(1), Truegd(2));

for aa = 1:length(estm)
    Basegd = xyz2gd(Base(aa,2:4));
    DDdXYZ(aa,:) = estm(aa,2:4) - Base(aa,2:4);
    DDdNEV(aa,:) = xyz2topo(DDdXYZ(aa,:), Basegd(1), Basegd(2));
    DDdN(aa,:) = DDdNEV(aa,1); 
    DDdE(aa,:) = DDdNEV(aa,2); 
    DDdV(aa,:) = DDdNEV(aa,3);
    DDdNE(aa,:) = (sqrt(DDdN(aa,1)^2 + DDdE(aa,1)^2));
    DDd3(aa,:) = (sqrt(DDdN(aa,1)^2 + DDdE(aa,1)^2 + DDdV(aa,1)^2));
    TrueDis(aa,1) = BsRvDis;
    TruedN(aa,1) = TruedNEV(1);
    TruedE(aa,1) = TruedNEV(2);
    TruedV(aa,1) = TruedNEV(3);
    
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
plot(DDdNE,'r-');
plot(DDd3,'b-');
plot(TrueDis,'k-');
drawnow
xlabel({['Double Differencing dNE = ', num2str(DDdNE_rms), '   std =', num2str(DDdNE_std)],...
        ['Double Differencing  3D = ', num2str(DDd3_rms), '   std =', num2str(DDd3_std)]}); 
ylabel('Distance(meter)');
% legend('2D(dNE) Distance','3D(d3D) Distance','True Distance')

subplot(3,2,1)
hold on; grid on;
xlim([0 length(estm)])
plot(DDdN,'r-');
plot(DDdE,'b-');
plot(DDdV,'k-');
xlabel('Epoch','fontsize',15);
ylabel('dN(meter)','fontsize',15);
legend('dN','dE','dV')

subplot(3,2,2)
hold on; grid on;
xlim([0 length(estm)])
plot(DDdN,'r-');
plot(TruedN,'k-');
xlabel({['Epoch'], ['Standard Deviation(dN) = ',num2str(DDdN_std)]});
ylabel('dN(meter)','fontsize',15);
legend('dN')

subplot(3,2,4)
hold on; grid on;
xlim([0 length(estm)])
plot(DDdE,'b-');
plot(TruedE,'k-');
xlabel({['Epoch'], ['Standard Deviation(dE) = ',num2str(DDdE_std)]});
ylabel('dE(meter)','fontsize',15);
legend('dE')

subplot(3,2,6)
hold on; grid on;
xlim([0 length(estm)])
plot(DDdV,'k-');
plot(TruedV,'k-');
xlabel({['Epoch'], ['Standard Deviation(dV) = ',num2str(DDdV_std)]});
ylabel('dV(meter)','fontsize',15);
legend('dV')