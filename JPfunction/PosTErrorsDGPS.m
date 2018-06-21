function [dXYZ, DGPSdXYZ, dNEV, DGPSdNEV] = PosTErrorsDgps(ttPlot, TruePos, EstPosT, EstPosTDGPS)
%
% function [dXYZ, DGPSdXYZ, dNEV, DGPSdNEV] = PosTErrorsDgps(tt, TruePos, EstPosT, EstPosT2)
%
% <input>   ttPlot: Time Tags [gs]
%           TruePos: True Position [1x3]
%           EstPosT: Estimated Positioni and Clock Error; c1-3(XYZ Estimates), c4(Clock Estimates)
%           EstPosT2: [DGPS]Estimated Positioni and Clock Error; c1-3(XYZ Estimates), c4(Clock Estimates)
%
% <output>  dXYZ: XYZ differences wrt TruePos
%           dNEV: NEV differences wrt TruePos
%           dXYZDGPS: DGPSXYZ differences wrt TruePos
%           dNEVDGPS: DGPSNEV differences wrt TruePos
%
%   Example : [dXYZ, DGPSdXYZ, dNEV, DGPSdNEV] = PosTErrorsDgps(tt, TruePos, EstPosT, EstPosT2)
%
%   Originally coded by Joonseong Gim, Jan 18, 2016
%% 입력받은 자료의 크기 및 참값에 해당하는 위경도 결정 결정
NoPos = length(EstPosT);
gd = xyz2gd(TruePos); TrueLat = gd(1); TrueLon = gd(2); 
%% 추정된 XYZ와 참값 XYZ의 차이값 결정
dXYZ = zeros(NoPos,3);
for k = 1:NoPos
    dXYZ(k,:) = EstPosT(k,1:3) - TruePos;
end
%% 추정된 DGPS XYZ와 참값 XYZ의 차이값 결정
DGPSdXYZ = zeros(NoPos,3);
for i = 1:NoPos
    DGPSdXYZ(i,:) = EstPosTDGPS(i,1:3) - TruePos;
end
%% dXYZ를 dNEV로 변환
dNEV = xyz2topo(dXYZ, TrueLat, TrueLon);
%% DGPSdXYZ를 DGPSdNEV로 변환
DGPSdNEV = xyz2topo(DGPSdXYZ, TrueLat, TrueLon);
%% 각 성분별 RMS 계산
dN = dNEV(:,1); dE = dNEV(:,2); dV = dNEV(:,3);
DGPSdN = DGPSdNEV(:,1); DGPSdE = DGPSdNEV(:,2); DGPSdV = DGPSdNEV(:,3); % DGPS
[dNE dU d3] = RMS(dNEV)
[DGPSdNE DGPSdU DGPSd3] = RMS(DGPSdNEV)

%% 그래프 그리는 작업 시작
rXY = max([abs(min(dNEV(:,1))) abs(max(dNEV(:,1))) abs(min(dNEV(:,2))) abs(max(dNEV(:,2)))]);
rXY = ceil(rXY);
figure(88);
subplot(3,4,[1,2,5,6])
hold on
plot(dE, dN,'o'); 
plot(DGPSdE, DGPSdN,'r.'); 
title('PP vs DGPS')
axis([-rXY rXY -rXY rXY]); grid on; 
axis square
xlabel({'\Delta E (meters)',...
    ['dNE = ', num2str(dNE),'  DGPSdNE = ', num2str(DGPSdNE)],...
    [' dV = ', num2str(dU),'   DGPSdV = ', num2str(DGPSdU)],...
    [' d3 = ', num2str(d3),'   DGPSd3 = ', num2str(DGPSd3)]}); 
ylabel('\Delta N (meters)')
%% 그래프 우측
tHour = mod(ttPlot, 86400);
tHour = tHour/3600;
dT = EstPosT(:,4);
DGPSdT = EstPosTDGPS(:,4);
subplot(3,4,[3,4])
plot(tHour, dN, '.r:', tHour, dE, '.b:'); axis([min(tHour) max(tHour) min(dN) max(dN)]); grid on;
legend('\Delta N', '\Delta E')
ylabel('\Delta H (meters)');
subplot(3,4,[7,8])
plot(tHour, DGPSdN, '.r:', tHour, DGPSdE, '.b:'); axis([min(tHour) max(tHour) min(DGPSdN) max(DGPSdN)]); grid on;
legend('DGPS \Delta N', 'DGPS \Delta E')
ylabel('\Delta H (meters)');
subplot(3,4,[11,12]);
plot(tHour, dV, '.:',tHour, DGPSdV, '.r:'); xlim([min(tHour) max(tHour)]); grid on;
legend('\Delta U', 'DGPS \Delta U')
ylabel('\Delta U (meters)')
% subplot(4,4,[15,16]);
% plot(tHour, dT, '.:', tHour, DGPSdT, '.r:'); axis([min(tHour) max(tHour) min(dT) max(dT)]); grid on; 
% ylabel('\delta t_r (meters)')
xlabel('Hours');



