function SPdopplot(QMfile)

% clear all
% close all
% load('pseudoranges_log_2017_09_14_17_57_01.mat');

qm = load(QMfile);
SVs = unique(qm(:,2));
CCC = 299792458;                                    % Speed of Light
WL = CCC/1575.42e6;


figure(1)
hold on; grid on;
xlim([min(qm(:,1)), max(qm(:,1))]);

for i=1:length(SVs)
    prn = SVs(i);
    GS = qm(find(qm(:,2) == prn & qm(:,3) == 120),1);
    PR = qm(find(qm(:,2) == prn & qm(:,3) == 120),:);
    PRR = qm(find(qm(:,2) == prn & qm(:,3) == 131),:);
    PRR(:,4) = PRR(:,4)*-WL;
    delGS = diff(GS);
    delPR = [GS(2:end), diff(PR(:,4))];
    y = delPR(:,2)./delGS;
    %     h = plot(delPR(:,1), delPR(:,2)); set(h, 'Marker','.','Markersize',4)
    h = plot(delPR(:,1), y); set(h, 'Marker','.','Markersize',4)
    ht=text(min(GS),y(1),int2str(prn),'Color',[.5 .5 .5]);
    set(ht,'HorizontalAlignment','right')
    colors = get(h,'Color');
    h2 = plot(PRR(:,1),PRR(:,4),'-k');
%     set(h2, 'Color',[.5 .5 .5])
    text(max(qm(:,1)),y(end),int2str(prn),'Color',colors);
    ts = ('diff(raw pr)/diff(time) and reported prr');
    title(ts)
    ylabel('(m/s)')
    xlabel('GPS second')
end

for i=1:length(SVs)
    figure(2)
    xlim([min(qm(:,1)), max(qm(:,1))]);
    subplot(length(SVs),1,i)
    hold on; grid on;
    xlim([min(qm(:,1)), max(qm(:,1))]);
    prn = SVs(i);
    GS = qm(find(qm(:,2) == prn & qm(:,3) == 120),1);
    PR = qm(find(qm(:,2) == prn & qm(:,3) == 120),:);
    PRR = qm(find(qm(:,2) == prn & qm(:,3) == 131),:);
    PRR(:,4) = PRR(:,4)*-WL;
    delGS = diff(GS);
    delPR = [GS(2:end), diff(PR(:,4))];
    y = delPR(:,2)./delGS;
    h = plot(delPR(:,1), y);set(h, 'Color',[.5 .5 .5])
    plot(PRR(:,1),PRR(:,4),'-b'); 
    text(max(qm(:,1)),PRR(end,4),int2str(prn),'Color',[.5 .5 .5]);
       
end

for i=1:length(SVs)
    figure(i+2)
    hold on; grid on;
    xlim([min(qm(:,1)), max(qm(:,1))]);
    prn = SVs(i);
    GS = qm(find(qm(:,2) == prn & qm(:,3) == 120),1);
    PR = qm(find(qm(:,2) == prn & qm(:,3) == 120),:);
    PRR = qm(find(qm(:,2) == prn & qm(:,3) == 131),:);
    PRR(:,4) = PRR(:,4)*-WL;
    delGS = diff(GS);
    delPR = [GS(2:end), diff(PR(:,4))];
    y = delPR(:,2)./delGS;
    h = plot(delPR(:,1), y);set(h, 'Color',[.5 .5 .5])
    plot(PRR(:,1),PRR(:,4),'-b');
    text(max(qm(:,1)),PRR(end,4),int2str(prn),'Color',[.5 .5 .5],'fontsize',20);
    ylabel('(m/s)','fontsize',20)
    xlabel('GPS second','fontsize',20)
    title(['PRN = ',num2str(prn)])
    legend('\Delta Pseudorange','PseudorangeRateMetersPerSecond')

    
end

