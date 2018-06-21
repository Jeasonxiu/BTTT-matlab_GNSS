function [dXYZ, dNEV] = PosTErrorsGPGGA(ttPlot, TruePos, EstPosT,name,yyyy,mm,dd,point)
%
% function [dXYZ, dNEV] = PosTErrorsGPGGA(tt, TruePos, EstPosT)
%
% <input>   ttPlot: Time Tags [gs]
%           TruePos: True Position [1x3]
%           EstPosT: Estimated Positioni and Clock Error; c1-3(XYZ Estimates), c4(Clock Estimates)
%
% <output>  dXYZ: XYZ differences wrt TruePos
%           dNEV: NEV differences wrt TruePos
%
% Copyright: Kwan-Dong Park, 9/20/2014
% Modified by JOON, 22/02/2016
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
%% �� ���к� RMS ���
dN = dNEV(:,1); dE = dNEV(:,2); dV = dNEV(:,3);
[dNE dU d3] = RMS(dNEV)
% dNE = sqrt(dN.^2 + dE.^2);        %rmsH = myRMS(dNE);
%                                   %rmsV = myRMS(dV);
% d3 = sqrt(dN.^2 + dE.^2 + dV.^2); %rms3 = myRMS(d3);
% fprintf('RMS Values \n H:%8.2f \n V:%8.2f \n d3:%7.2f\n', rmsH, rmsV, rms3)
%% �׷��� �׸��� �۾� ����
rXY = max([abs(min(dNEV(:,1))) abs(max(dNEV(:,1))) abs(min(dNEV(:,2))) abs(max(dNEV(:,2)))]);
rXY = ceil(rXY);
figure(66)
subplot(2,4,[1,2,5,6])
plot(dE, dN,'o'); 
axis([-rXY rXY -rXY rXY]); grid on; 
axis square
% title({[name,'(point: ',point,')' ]; [num2str(yyyy),'-',num2str(mm),'-',num2str(dd)];...
%     [' epoch(', num2str(length(EstPosT)),')']});
xlabel({'\Delta E (meters)',...
    ['dNE = ', num2str(dNE)],...
    [' dV = ', num2str(dU)],...
    [' 3D = ', num2str(d3)]}); ; ylabel('\Delta N (meters)')
%% �׷��� ����
tHour = mod(ttPlot, 86400);
tHour = tHour/3600;
    subplot(2,4,[3,4])
    plot(tHour, dN, '.r:', tHour, dE, '.b:'); 
    xlim([min(tHour) max(tHour)]);
%     axis([min(tHour) max(tHour) min(dN) max(dN)]); 
    grid on;
    legend('\Delta N', '\Delta E')
    ylabel('\Delta H (meters)');
    subplot(2,4,[7,8])
    plot(tHour, dV, '.:'); axis([min(tHour) max(tHour) min(dV) max(dV)]); grid on;
    ylabel('\Delta U (meters)')
    xlabel('Hours');


