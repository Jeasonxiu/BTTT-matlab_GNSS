function [dXYZ, dNEV] = PosTErrorsNMEA(nmea, TruePos)
%
% [dXYZ, dNEV] = PosTErrorsNMEA(nmea, TruePos)
%
% <input>   nmea: nmez array : after WriteNMEA
%           TruePos: True Position [1x3](ECEF)

%
% <output>  dXYZ: XYZ differences wrt TruePos
%           dNEV: NEV differences wrt TruePos
%
% Copyright: Joonseong Gim, 3/29/2018
%
%% 입력받은 자료의 크기 및 참값에 해당하는 위경도 결정 결정
FinalTTs = unique(nmea(:,1));
%% 입력받은 자료의 크기 및 참값에 해당하는 위경도 결정 결정
NoPos = length(FinalTTs);
gd = xyz2gd(TruePos); TrueLat = gd(1); TrueLon = gd(2); 
%% 추정된 XYZ와 참값 XYZ의 차이값 결정
dXYZ = zeros(NoPos,3);
for k = 1:NoPos
    gs = FinalTTs(k);
    gd = nmea(find(nmea(:,1) == gs),2:4);
    gd = unique(gd);
    xyz = gd2xyz(gd);
    dXYZ(k,:) = xyz - TruePos;
end
%% dXYZ를 dNEV로 변환
dNEV = xyz2topo(dXYZ, TrueLat, TrueLon);
absdNEV = abs(dNEV); Max_error = max(absdNEV(:,1:2)); Max_error = max(Max_error);
%% 각 성분별 RMS 계산
dN = dNEV(:,1); dE = dNEV(:,2); dV = dNEV(:,3);
[dNE dU d3] = RMS(dNEV);
% dNE = sqrt(dN.^2 + dE.^2);        %rmsH = myRMS(dNE);
%                                   %rmsV = myRMS(dV);
% d3 = sqrt(dN.^2 + dE.^2 + dV.^2); %rms3 = myRMS(d3);
% fprintf('RMS Values \n H:%8.2f \n V:%8.2f \n d3:%7.2f\n', rmsH, rmsV, rms3)
%% 그래프 그리는 작업 시작
rXY = max([abs(min(dNEV(:,1))) abs(max(dNEV(:,1))) abs(min(dNEV(:,2))) abs(max(dNEV(:,2)))]);
rXY = ceil(rXY);
figure();
subplot(2,4,[1,2,5,6])
plot([-200,200],[0,0],'r-');hold on; grid on;
plot([0,0],[-200,200],'r-')
plot(dE, dN,'bo')

% axis([-rXY rXY -rXY rXY]);
% axis([-5 5 -5 5]);
axis([-Max_error Max_error -Max_error Max_error]);
% axis([-ceil(dNE) ceil(dNE) -ceil(dNE) ceil(dNE)]);
axis square
grid on;
xlabel({'\Delta E (meters)',...
    ['dNE = ', num2str(decimal(dNE))],...
    [' dV = ', num2str(decimal(dU))],...
    [' 3D = ', num2str(decimal(d3))]}); ; ylabel('\Delta N (meters)')
%% 그래프 우측
tHour = mod(FinalTTs, 86400);
tHour = tHour/3600;
subplot(2,4,[3,4])
plot(tHour, dN, '.r-', tHour, dE, '.b-'); axis([min(tHour) max(tHour) min(dN) max(dN)]); grid on;
ylim([-Max_error Max_error]);
legend('\Delta N', '\Delta E')
ylabel('\Delta N, \Delta E (meters)');
subplot(2,4,[7,8])
plot(tHour, dV, 'b.-'); axis([min(tHour) max(tHour) min(dV) max(dV)]); grid on;
% ylim([-Max_error Max_error]);
ylabel('\Delta V (meters)')
fprintf('H RMSE = %3.2f cm : V RMSE = %3.2f cm : 3D RMSE = %3.2f cm \n',dNE, dU, d3);
