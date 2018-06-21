% function [dXYZ, dNEV] = PosTErrorsGPGGA2(ttPlot, TruePos, EstPosT)
%
% function [dXYZ, dNEV] = PosTErrorsGPGGA2(tt, TruePos, EstPosT)
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
%% 입력받은 자료의 크기 및 참값에 해당하는 위경도 결정 결정
NoPos = length(EstPosT);
gd = xyz2gd(TruePos); TrueLat = gd(1); TrueLon = gd(2); 
%% 추정된 XYZ와 참값 XYZ의 차이값 결정
dXYZ = zeros(NoPos,3);
for k = 1:NoPos
    dXYZ(k,:) = EstPosT(k,1:3) - TruePos;
    
end
% dXYZ(:,4) = EstPosT(:,4);
%% dXYZ를 dNEV로 변환
dNEV = xyz2topo(dXYZ, TrueLat, TrueLon);
%% 각 성분별 RMS 계산
% dN = dNEV(:,1); dE = dNEV(:,2); dV = dNEV(:,3);
dN(:,1) = dNEV(:,1); dE(:,1) = dNEV(:,2); dV(:,1) = dNEV(:,3);
dN(:,2) = EstPosT(:,4); dE(:,2) = EstPosT(:,4); dV(:,2) = EstPosT(:,4);
GGA = [dN(find(dN(:,2) == 1)), dE(find(dE(:,2) == 1))];
FLOAT = [dN(find(dN(:,2) == 5)), dE(find(dE(:,2) == 5))];
FIXED = [dN(find(dN(:,2) == 4)), dE(find(dE(:,2) == 4))];
[dNE dU d3] = RMS(dNEV)
% dNE = sqrt(dN.^2 + dE.^2);        %rmsH = myRMS(dNE);
%                                   %rmsV = myRMS(dV);
% d3 = sqrt(dN.^2 + dE.^2 + dV.^2); %rms3 = myRMS(d3);
% fprintf('RMS Values \n H:%8.2f \n V:%8.2f \n d3:%7.2f\n', rmsH, rmsV, rms3)
%% 그래프 그리는 작업 시작
rXY = max([abs(min(dNEV(:,1))) abs(max(dNEV(:,1))) abs(min(dNEV(:,2))) abs(max(dNEV(:,2)))]);
rXY = ceil(rXY);
figure(66)
subplot(2,4,[1,2,5,6])
hold on
plot(GGA(:,2), GGA(:,1),'bo'); 
plot(FLOAT(:,2), FLOAT(:,1),'go'); 
plot(FIXED(:,2), FIXED(:,1),'ro'); 
axis([-rXY rXY -rXY rXY]); grid on; 
axis square
% title({[name,'(point: ',point,')' ]; [num2str(yyyy),'-',num2str(mm),'-',num2str(dd)];...
%     [' epoch(', num2str(length(EstPosT)),')']});
xlabel({'\Delta E (meters)',...
    ['dNE = ', num2str(dNE)],...
    [' dV = ', num2str(dU)],...
    [' 3D = ', num2str(d3)]}); ; ylabel('\Delta N (meters)')
%% 그래프 우측
tHour = mod(ttPlot, 86400);
tHour = tHour/3600;
    subplot(2,4,[3,4])
    plot(tHour, dN(:,1), '.r:', tHour, dE(:,1), '.b:'); 
    xlim([min(tHour) max(tHour)]);
%     axis([min(tHour) max(tHour) min(dN) max(dN)]); 
    grid on;
    legend('\Delta N', '\Delta E')
    ylabel('\Delta H (meters)');
    subplot(2,4,[7,8])
    plot(tHour, dV(:,1), '.:'); axis([min(tHour) max(tHour) min(dV(:,1)) max(dV(:,1))]); grid on;
    ylabel('\Delta U (meters)')
    xlabel('Hours');


