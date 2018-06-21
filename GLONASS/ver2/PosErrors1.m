function [dXYZ, dNEV] = PosErrors1(ttPlot, TruePos, EstPos,visiSat)

% �⺻�� PosError�� ����, ���������� �߰�
% Mi-So Kim, April 13, 2015 

%% �Է¹��� �ڷ��� ũ�� �� ������ �ش��ϴ� ���浵 ���� ����
NoPos = length(EstPos);
gd = xyz2gd(TruePos); TrueLat = gd(1); TrueLon = gd(2); 

%% ������ XYZ�� ���� XYZ�� ���̰� ����
dXYZ = zeros(NoPos,3);
for k = 1:NoPos
    dXYZ(k,:) = EstPos(k,:) - TruePos;
end

%% dXYZ�� dNEV�� ��ȯ
dNEV = xyz2topo_all(dXYZ, TrueLat, TrueLon); % ��� dXYZ�� xyz2topo�� �̿��Ͽ� dNEV�� ��ȯ�ߴµ� ����� �ƴϰ� �ϳ��� ������ ����.

%% �� ���к� RMS ���
% dN = dNEV(:,1); dE = dNEV(:,2); dV = dNEV(:,3);
% dNE = sqrt(dN.^2 + dE.^2);        rmsH = std(dNE);
%                                   rmsV = std(dV);
% d3 = sqrt(dN.^2 + dE.^2 + dV.^2); rms3 = std(d3);

dN = dNEV(:,1); dE = dNEV(:,2); dV = dNEV(:,3);
dNE = sqrt(dN.^2 + dE.^2);        rmsH = myRMS(dNE);
                                  rmsV = myRMS(dV);
d3 = sqrt(dN.^2 + dE.^2 + dV.^2); rms3 = myRMS(d3);
fprintf('RMS Values \n H:%8.2f \n V:%8.2f \n 3D:%7.2f\n', rmsH, rmsV, rms3)

%% �׷��� �׸��� �۾� ����
rXY = max([abs(min(dNEV(:,1))) abs(max(dNEV(:,1))) abs(min(dNEV(:,2))) abs(max(dNEV(:,2)))]);
rXY = ceil(rXY);
subplot(3,3,4)
plot(dE, dN,'o'); 
axis([-15 15 -15 15]); grid on; 
xlabel('\Delta E (meters)'); ylabel('\Delta N (meters)')

tHour = mod(ttPlot, 86400);
tHour = tHour/3600;

subplot(3,3,2:3) % N & E
plot(tHour, dN, '.:b',tHour,dE,'.:r'); axis([0 max(tHour) -10 10]); grid on;
ylabel('\Delta N E (meters)'); legend('dN','dE');

subplot(3,3,5:6); % V
plot(tHour, dV, '.:'); axis([0 max(tHour) min(dV) max(dV)]); grid on; 
ylabel('\Delta V (meters)')

subplot(3,3,8:9); % ��������
stairs(tHour, visiSat(:,1), '-.ok'); axis tight; grid on;
xlabel('Hours'); ylabel('Visiablity Satellite');

% subplot(3,3,2:3)
% plot(tHour, dN, '.:'); axis([0 max(tHour) min(dN) max(dN)]); grid on;
% ylabel('\Delta N (meters)')
% 
% subplot(3,3,5:6);
% plot(tHour, dE, '.:'); axis([0 max(tHour) min(dE) max(dE)]); grid on; 
% ylabel('\Delta E (meters)')
% 
% subplot(3,3,8:9);
% plot(tHour, dV, '.:'); axis([0 max(tHour) min(dV) max(dV)]); grid on;
% xlabel('Hours'); ylabel('\Delta U (meters)')

%% �߰��۾�
%: �ΰ��� �۾�. �߽ɼ� �׸���, ���Ϳ��� ��Į�� ����(�ڸǵ忡�� ��)
