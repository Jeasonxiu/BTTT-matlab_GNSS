function [dXYZ, dNEV] = PosTErrorsCP(ttPlot, TruePos, EstPosT, Estcp)
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
%
%% �Է¹��� �ڷ��� ũ�� �� ������ �ش��ϴ� ���浵 ���� ����
NoPos = length(EstPosT);
gd = xyz2gd(TruePos); TrueLat = gd(1); TrueLon = gd(2); 
%% ������ XYZ�� ���� XYZ�� ���̰� ����
dXYZ = zeros(NoPos,3);
dXYZcp = zeros(NoPos,3);
for k = 1:NoPos
    dXYZ(k,:) = EstPosT(k,1:3) - TruePos;
    dXYZcp(k,:) = Estcp(k,1:3) - TruePos;
end
%% dXYZ�� dNEV�� ��ȯ
dNEV = xyz2topo(dXYZ, TrueLat, TrueLon);
dNEVcp = xyz2topo(dXYZcp, TrueLat, TrueLon);
% absdNEV = abs(dNEV); Max_error = max(absdNEV(:,1:2)); Max_error = max(Max_error);
%% �� ���к� RMS ���
dN = dNEV(:,1); dE = dNEV(:,2); dV = dNEV(:,3);
dNcp = dNEVcp(:,1); dEcp = dNEVcp(:,2); dVcp = dNEVcp(:,3);
[dNE dU d3] = RMS(dNEV)
[dNEcp dUcp d3cp] = RMS(dNEVcp)
% dNE = sqrt(dN.^2 + dE.^2);        %rmsH = myRMS(dNE);
%                                   %rmsV = myRMS(dV);
% d3 = sqrt(dN.^2 + dE.^2 + dV.^2); %rms3 = myRMS(d3);
% fprintf('RMS Values \n H:%8.2f \n V:%8.2f \n 3D:%7.2f\n', rmsH, rmsV, rms3)
%% �׷��� �׸��� �۾� ����
% rXY = max([abs(min(dNEV(:,1))) abs(max(dNEV(:,1))) abs(min(dNEV(:,2))) abs(max(dNEV(:,2)))]);
% rXY = ceil(rXY);
figure();
subplot(1,2,1)
plot(dE, dN,'go', dEcp, dNcp,'r+'); 
% axis([-rXY rXY -rXY rXY]); 
% axis([-5 5 -5 5]);
axis([-10 10 -10 10]);
% axis([-ceil(dNE) ceil(dNE) -ceil(dNE) ceil(dNE)]);
axis square
grid on; 
legend('Smartphone', 'DGPS-CP')
%% �׷��� ����
tHour = mod(ttPlot, 86400);
tHour = tHour/3600;
dT = EstPosT(:,4);
subplot(1,2,2)
plot(tHour, dV, 'go:', tHour, dVcp, 'r+:'); 
ylim([-20 20]);
axis square
grid on;
legend('Smartphone', 'DGPS-CP')
