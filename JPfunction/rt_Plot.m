function rt_Plot(FinalTTs, gs, EstPos, TruePos)
%
% function rt_Plot(FinalTTs, gs, estm, TruePos)
%
% DO: Plot real-time positional error 
%
% <input>   FinalTTs : 전체 수집데이터 시간           
%           gs
%           EstPos: Estimated position in XYZ
%           TruePos: True position in XYZ
%
% <output>  Real-time Plot
%
% Copyright: Joonseong Gim, 4/12/2018 @PPSoln

% FinalTTs to hour
tHour = mod(FinalTTs, 86400); tHour = tHour/3600;
% TruePos geodetic
TruePos_gd = xyz2gd(TruePos);

% first epoch base plot
if gs == FinalTTs(1)
    % real-time plot animation parameter
    display = get(0,'MonitorPositions'); 
    figure(111)
    set(gcf,'position',[display(1,3)-600, display(1,4)-651, 580, 631]);
    subplot(5,1,[1:4]);
    polarplot(0,0,'r+','Markersize',10); hold on; grid on;
    P_Labels = {'E';[];[];'N';[];[];'W';[];[];'S';[];[]};
    set(gca, 'ThetaTickLabel', P_Labels);
    set(gca, 'RAxisLocation', 0);
    rlim([0, 0.5]);
    set(gca, 'rticklabels', {'0';'0.1';'0.2';'0.3';'0.4';'0.5 m'});
    subplot(5,1,5);
    hold on; grid on;
    xlim([min(tHour) max(tHour)])
    drawnow limitrate
end

% gs to hour
gs_hour = mod(gs, 86400); gs_hour = gs_hour/3600;
%% dXYZ, dNEV
dXYZ = EstPos - TruePos;
dNEV = [gs, xyz2topo(dXYZ, TruePos_gd(1), TruePos_gd(2))];
% real-time plot
figure(111)
subplot(5,1,[1:4])
title(['H error : ',num2str(norm(dNEV(2:3)),'%4.2f'),'(m)'])
polarplot(atan2(dNEV(2),dNEV(3)), norm(dNEV(2:3)),'bo','Markersize',5,...
    'markerfacecolor','r');
rlimmax = 0.5;
if norm(dNEV(2:3)) > rlimmax
    rlimmax = ceil(norm(dNEV(2:3))*10)/10;
    rlim([0, rlimmax]);
    Rtick = get(gca, 'rtick');
    set(gca, 'rticklabels', Rtick);
    Rtick = get(gca, 'rticklabels');
    Rtick = strcat(Rtick,' m');
    set(gca, 'rticklabels', Rtick);
end

% if norm(dNEV(2:3)) < 0.5
%     rlim([0, 0.5]);
%     set(gca, 'rticklabels', {'0';'0.1';'0.2';'0.3';'0.4';'0.5 m'});
%     polarplot(atan2(dNEV(2),dNEV(3)), norm(dNEV(2:3)),'bo','Markersize',5,...
%     'markerfacecolor','r');
% elseif norm(dNEV(2:3)) > 0.5 &  norm(dNEV(2:3)) < 1
%     rlim([0, 1]);
%     set(gca, 'rticklabels', {'0';'0.2';'0.4';'0.6';'0.8';'1 m'});
%     polarplot(atan2(dNEV(2),dNEV(3)), norm(dNEV(2:3)),'bo','Markersize',5,...
%     'markerfacecolor','y');
% elseif norm(dNEV(2:3)) > 1 &  norm(dNEV(2:3)) < 2
%     rlim([0, 1]);
%     set(gca, 'rticklabels', {'0';'0.4';'0.8';'1.2';'1.6';'2 m'});
%     polarplot(atan2(dNEV(2),dNEV(3)), norm(dNEV(2:3)),'bo','Markersize',5,...
%     'markerfacecolor','g');
% elseif norm(dNEV(2:3)) > 2
%     rlim([0, ceil(norm(dNEV(2:3)))]);
%     Rtick = get(gca, 'rticklabels');
%     Rtick{end} = strcat(Rtick{end},' m');
%     set(gca, 'rticklabels', Rtick);
%     polarplot(atan2(dNEV(2),dNEV(3)), norm(dNEV(2:3)),'bo','Markersize',5,...
%     'markerfacecolor','k');
% end
subplot(5,1,5)
title(['V error : ',num2str(dNEV(4),'%4.2f'),'(m)'])
plot(gs_hour, dNEV(4), 'bo','Markersize',5,'markerfacecolor','r')
xlabel('Time Elapsed (Hours)');
