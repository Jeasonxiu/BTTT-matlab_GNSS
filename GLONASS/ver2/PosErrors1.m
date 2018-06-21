function [dXYZ, dNEV] = PosErrors1(ttPlot, TruePos, EstPos,visiSat)

% 기본의 PosError을 수정, 가시위성을 추가
% Mi-So Kim, April 13, 2015 

%% 입력받은 자료의 크기 및 참값에 해당하는 위경도 결정 결정
NoPos = length(EstPos);
gd = xyz2gd(TruePos); TrueLat = gd(1); TrueLon = gd(2); 

%% 추정된 XYZ와 참값 XYZ의 차이값 결정
dXYZ = zeros(NoPos,3);
for k = 1:NoPos
    dXYZ(k,:) = EstPos(k,:) - TruePos;
end

%% dXYZ를 dNEV로 변환
dNEV = xyz2topo_all(dXYZ, TrueLat, TrueLon); % 행렬 dXYZ을 xyz2topo를 이용하여 dNEV로 변환했는데 행렬이 아니고 하나의 값으로 나옴.

%% 각 성분별 RMS 계산
% dN = dNEV(:,1); dE = dNEV(:,2); dV = dNEV(:,3);
% dNE = sqrt(dN.^2 + dE.^2);        rmsH = std(dNE);
%                                   rmsV = std(dV);
% d3 = sqrt(dN.^2 + dE.^2 + dV.^2); rms3 = std(d3);

dN = dNEV(:,1); dE = dNEV(:,2); dV = dNEV(:,3);
dNE = sqrt(dN.^2 + dE.^2);        rmsH = myRMS(dNE);
                                  rmsV = myRMS(dV);
d3 = sqrt(dN.^2 + dE.^2 + dV.^2); rms3 = myRMS(d3);
fprintf('RMS Values \n H:%8.2f \n V:%8.2f \n 3D:%7.2f\n', rmsH, rmsV, rms3)

%% 그래프 그리는 작업 시작
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

subplot(3,3,8:9); % 가시위성
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

%% 추가작업
%: 두가지 작업. 중심선 그리기, 벡터에서 스칼라 빼기(코맨드에선 됨)
