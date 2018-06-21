function [dXYZ, dNEV] = PosTErrors2(ttPlot, TruePos, EstPosT,Sats)
%
% function [dXYZ, dNEV] = PosTErrors(tt, TruePos, EstPosT)
%
% <input>   ttPlot: Time Tags [gs]
%           TruePos: True Position [1x3]
%           EstPosT: Estimated Positioni and Clock Error; c1-3(XYZ Estimates), c4(Clock Estimates)
%
% <output>  dXYZ: XYZ differences wrt TruePos
%           dNEV: NEV differences wrt TruePos
%
% Copyright: Kwan-Dong Park, 9/20/2014
% 교수님 코드를 기반으로 RMSE가 같이 표시되도록 수정

%% 입력받은 자료의 크기 및 참값에 해당하는 위경도 결정 결정
NoPos = length(EstPosT);
gd = xyz2gd(TruePos); TrueLat = gd(1); TrueLon = gd(2); 
%% 추정된 XYZ와 참값 XYZ의 차이값 결정
dXYZ = zeros(NoPos,3);
for k = 1:NoPos
    dXYZ(k,:) = EstPosT(k,1:3) - TruePos;
end
%% dXYZ를 dNEV로 변환
dNEV = xyz2topo(dXYZ, TrueLat, TrueLon);
% absdNEV = abs(dNEV); Max_error = max(absdNEV(:,1:2)); Max_error = max(Max_error);
%% 각 성분별 RMS 계산
dN = dNEV(:,1); dE = dNEV(:,2); dV = dNEV(:,3);
dNE = sqrt(dN.^2 + dE.^2);        %rmsH = myRMS(dNE);
dU = sqrt(dV.^2);                 %rmsV = myRMS(dV);
d3 = sqrt(dN.^2 + dE.^2 + dV.^2); %rms3 = myRMS(d3);
% fprintf('RMS Values \n H:%8.2f \n V:%8.2f \n 3D:%7.2f\n', rmsH, rmsV, rms3)
%% 그래프 그리는 작업 시작
rXY = max([abs(min(dNEV(:,1))) abs(max(dNEV(:,1))) abs(min(dNEV(:,2))) abs(max(dNEV(:,2)))]);
rXY = ceil(rXY);
figure();
subplot(3,4,[1,2,5,6])
plot(dE, dN,'o'); 
axis([-rXY rXY -rXY rXY]); 
% axis([-5 5 -5 5]);
% axis([-Max_error Max_error -Max_error Max_error]);
% axis([-ceil(dNE) ceil(dNE) -ceil(dNE) ceil(dNE)]);
axis square
grid on; 
xlabel({'\Delta E (meters)',...
    ['dNE = ', num2str(rms(dNE)),...
    '    dV = ', num2str(rms(dU)),...
    '    3D = ', num2str(rms(d3))]}); ; ylabel('\Delta N (meters)')
%% 그래프 우측
tHour = mod(ttPlot, 86400);
tHour = tHour/3600;
dT = EstPosT(:,4);
subplot(3,4,[3,4])
plot(tHour, dN, '.r:', tHour, dE, '.b:'); 
% axis([min(tHour) max(tHour) min(dN) max(dN)]); grid on;
xlim([min(tHour) max(tHour)]); grid on;
legend('\Delta N', '\Delta E')
ylabel('\Delta H (meters)');
subplot(3,4,[7,8])
plot(tHour, dV, '.:'); 
% axis([min(tHour) max(tHour) min(dV) max(dV)]); grid on;
xlim([min(tHour) max(tHour)]); grid on;
ylabel('\Delta V (meters)')
subplot(3,4,[11,12])
plot(tHour, dT, '.:'); 
% axis([min(tHour) max(tHour) min(dT) max(dT)]); grid on; 
xlim([min(tHour) max(tHour)]); grid on;
ylabel('\delta t_r (meters)')
xlabel('Hours');

if nargin == 4
    subplot(3,4,[9,10])
    plot(tHour, Sats(:,1), 'r-'); axis([min(tHour) max(tHour) 4 max(Sats(:,1))]); 
    hold on; grid on;
    plot(tHour, Sats(:,2), 'b-'); axis([min(tHour) max(tHour) 4 max(Sats(:,1))]); 
    ylabel('Number of Sats')
    xlabel('Hours');
    legend('Obsevation','Used')
end
disp(    ['dNE = ', num2str(rms(dNE)),...
    '    dV = ', num2str(rms(dU)),...
    '    3D = ', num2str(rms(d3))])