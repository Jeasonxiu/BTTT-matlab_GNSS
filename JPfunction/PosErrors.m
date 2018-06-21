function [dXYZ, dNEV] = PosErrors(ttPlot, TruePos, EstPos, Scount)

% load EstPos.mat
% 
% ttPlot = EstPos(:,1);
% EstPos = EstPos(:,2:4);

% TruePos = [-3.0267    4.0672    3.8572].*1.0e+006;

%% 입력받은 자료의 크기 및 참값에 해당하는 위경도 결정 결정
NoPos = length(EstPos);
gd = xyz2gd(TruePos); TrueLat = gd(1); TrueLon = gd(2); 

%% 추정된 XYZ와 참값 XYZ의 차이값 결정
dXYZ = zeros(NoPos,3);
for k = 1:NoPos
    dXYZ(k,:) = EstPos(k,1:3) - TruePos;
end

%% dXYZ를 dNEV로 변환
dNEV = xyz2topo(dXYZ, TrueLat, TrueLon); % 행렬 dXYZ을 xyz2topo를 이용하여 dNEV로 변환했는데 행렬이 아니고 하나의 값으로 나옴.
% dNEV = xyz2topo(dXYZ, TrueLat, TrueLon);
%% 각 성분별 RMS 계산
dN = dNEV(:,1); dE = dNEV(:,2); dV = dNEV(:,3);
% dNE = sqrt(dN.^2 + dE.^2);        rmsH = std(dNE);
%                                   rmsV = std(dV);
% d3 = sqrt(dN.^2 + dE.^2 + dV.^2); rms3 = std(d3);

dN = dNEV(:,1); dE = dNEV(:,2); dV = dNEV(:,3);
dNE = sqrt(dN.^2 + dE.^2);        rmsH = sqrt(mean(dNE.^2));
                                  rmsV = sqrt(mean(dV.^2));
d3 = sqrt(dN.^2 + dE.^2 + dV.^2); rms3 = sqrt(mean(d3.^2));
fprintf('RMS Values \n H:%8.2f \n V:%8.2f \n 3D:%7.2f\n', rmsH, rmsV, rms3)

%% 그래프 그리는 작업 시작
figure(101)
rXY = max([abs(min(dNEV(:,1))) abs(max(dNEV(:,1))) abs(min(dNEV(:,2))) abs(max(dNEV(:,2)))]);
rXY = ceil(rXY);
subplot(3,3,4)
plot(dE, dN,'o:'); 
axis([-rXY rXY -rXY rXY]); grid on; 
% axis([-30 30 -30 30]); 
grid on; 
xlabel('\Delta E (meters)'); ylabel('\Delta N (meters)')

tHour = mod(ttPlot, 86400);
tHour = tHour/3600;

subplot(3,3,2:3)
plot(tHour, dN, '.:b'); hold on; 
plot(tHour, dE, '.:r'); hold on;
legend('\Delta N', '\Delta E')
% plot(tHour, dNE, '.:g'); 
axis([min(tHour) max(tHour) min([dN;dE]) max([dN;dE])]);
% axis([min(tHour) max(tHour) -10 10]); 
grid on;
ylabel('\Delta N E(meters)')

% subplot(3,3,5:6);
% plot(tHour, dE, '.:'); axis([min(tHour) max(tHour) min(dE) max(dE)]); grid on; 
% ylabel('\Delta E (meters)')

subplot(3,3,5:6);
plot(tHour, dV, '.:'); axis([min(tHour) max(tHour) min(dV) max(dV)]); grid on;
% plot(tHour, dV, '.:'); 
% axis([min(tHour) max(tHour) -20 20]); 
% grid on;
% xlabel('Hours'); 
ylabel('\Delta U (meters)')

subplot(3,3,8:9);
% plot(tHour, dV, '.:'); axis([min(tHour) max(tHour) min(dV) max(dV)]); grid on;
plot(tHour, Scount, '.:'); 
axis([min(tHour) max(tHour) min(Scount)-3 max(Scount)+3]); 
grid on;
xlabel('Hours');  ylabel('# of Sats');% ylabel('\Delta U (meters)')
