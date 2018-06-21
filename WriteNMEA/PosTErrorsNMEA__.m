function [dXYZ, dNEV] = PosTErrorsNMEA(ttPlot, TruePos, EstPosT, UsedSats)
%
% [dXYZ, dNEV] = PosTErrorsNMEA(ttPlot, TruePos, EstPosT, UsedSats)
%
% <input>   ttPlot: Time Tags [gs]
%           TruePos: True Position [1x3]
%           EstPosT: Estimated Position
%           UsedSats : Number of Sats for Positioing
%
% <output>  dXYZ: XYZ differences wrt TruePos
%           dNEV: NEV differences wrt TruePos
%
% Copyright: Kwan-Dong Park, 9/20/2014
% Modified by JOON, 02/01/2018
%
%% �Է¹��� �ڷ��� ũ�� �� ������ �ش��ϴ� ���浵 ���� ����
NoPos = length(EstPosT);
gd = xyz2gd(TruePos); TrueLat = gd(1); TrueLon = gd(2); 
%% ������ XYZ�� ���� XYZ�� ���̰� ����
dXYZ = zeros(NoPos,3);
for k = 1:NoPos
    dXYZ(k,:) = EstPosT(k,1:3) - TruePos;
end
%% dXYZ�� dNEV�� ��ȯ
dNEV = xyz2topo(dXYZ, TrueLat, TrueLon);
absdNEV = abs(dNEV); Max_error = max(absdNEV(:,1:2)); Max_error = max(Max_error);
%% �� ���к� RMS ���
dN = dNEV(:,1); dE = dNEV(:,2); dV = dNEV(:,3);
[dNE dU d3] = RMS(dNEV);
% dNE = sqrt(dN.^2 + dE.^2);        %rmsH = myRMS(dNE);
%                                   %rmsV = myRMS(dV);
% d3 = sqrt(dN.^2 + dE.^2 + dV.^2); %rms3 = myRMS(d3);
% fprintf('RMS Values \n H:%8.2f \n V:%8.2f \n d3:%7.2f\n', rmsH, rmsV, rms3)
%% �׷��� �׸��� �۾� ����
rXY = max([abs(min(dNEV(:,1))) abs(max(dNEV(:,1))) abs(min(dNEV(:,2))) abs(max(dNEV(:,2)))]);
rXY = ceil(rXY);
figure();
subplot(4,4,[1,2,5,6])
plot([-200,200],[0,0],'r-');hold on; grid on;
plot([0,0],[-200,200],'r-')
plot(dE, dN,'bo')

axis([-Max_error Max_error -Max_error Max_error]);
axis square
grid on;
xlabel({'\Delta E (meters)',...
    ['dNE = ', num2str(decimal(dNE))],...
    [' dV = ', num2str(decimal(dU))],...
    [' 3D = ', num2str(decimal(d3))]}); ; ylabel('\Delta N (meters)')
%% �׷��� ����
tHour = mod(ttPlot, 86400);
tHour = tHour/3600;
subplot(4,4,[3,4])
plot(tHour, dN, '.r-', tHour, dE, '.b-'); axis([min(tHour) max(tHour) min(dN) max(dN)]); grid on;
ylim([-Max_error Max_error]);
legend('\Delta N', '\Delta E')
ylabel('\Delta N, \Delta E (meters)');
subplot(4,4,[7,8])
plot(tHour, dV, 'b.-'); axis([min(tHour) max(tHour) min(dV) max(dV)]); grid on;
% ylim([-Max_error Max_error]);
ylabel('\Delta V (meters)')

if ~isempty(find(UsedSats(:,:) ~= 0))
    subplot(4,4,[11,12])
    stairs(tHour, UsedSats(:,1), '-ob');
    xlim([min(tHour) max(tHour)]); grid on;
    ylabel('Used Sats(Total)')
    subplot(4,4,[15,16])
    stairs(tHour, UsedSats(:,2), 'b.-'); grid on; hold on;
    stairs(tHour, UsedSats(:,3), 'r.-');
    stairs(tHour, UsedSats(:,4), 'c.-');
    stairs(tHour, UsedSats(:,5), 'k.-');
    xlim([min(tHour) max(tHour)]); grid on;
    ylabel('Used Sats')
    legend('GPS','BDS','GLO','QZSS','Location','northoutside','Orientation','horizontal')
end
fprintf('H RMSE = %3.2f cm : V RMSE = %3.2f cm : 3D RMSE = %3.2f cm \n',dNE, dU, d3);
