% function Skyplotnmea(nmeafile)
nmeafile =  'jprA214f_s7.txt';
[NMEAQM,nmeaqm] = writeNMEA2(nmeafile);

%% Plot Sky
% Az = NMEAQM(:,3);
% El = NMEAQM(:,4);
cutoff = 5;


figure(100);
hold on;
Az = NMEAQM(:,7); El = NMEAQM(:,6); Sg = NMEAQM(:,8);
SgSize = 40;
SgColor = Sg;


%    Draw Circle
    circlea = 0:pi/30:2*pi;
    cx = cos(circlea);
    cy = sin(circlea);
 
    for i= [30 60 90]
        plot(cx*i, cy*i, '-', 'color', 'k', 'linewidth', 1);
    end
 
    for i=[15 45 75]
        plot(cx*i, cy*i, ':', 'color', 'k', 'linewidth', 1);
    end

%    Draw Lines inside a Circle
    lenmax = 90;
    circleTick = (1:6)*pi/6;
    cosct = cos(circleTick); 
    sinct = sin(circleTick);
    cax = [-cosct; cosct];
    say = [-sinct; sinct];
    plot(lenmax*cax, lenmax*say, '-', 'color', 'k', 'linewidth', 1);

    xx = (El-90) .* -(sin(Az*pi/180));
    yy = (El-90) .* -(cos(Az*pi/180));

%    Draw point Signal Strength value
    scatter(xx, yy, SgSize, SgColor,'filled');
    title('SkyPlot');
    colormap(jet);
    colorbar('location', 'EastOutside');
    caxis([0 60]);
    
%    Insert direction text
    rlen = 1.06 * lenmax;
    for i = 1 : length(circleTick) 
        ticm1 = int2str(i*30);
        ticm2 = int2str(180+i*30);
        if ticm2 == '360'
            ticm2 =' ';
            %ticm2 ='N';
        end
        text( rlen*sinct(i),  rlen * cosct(i), ticm1, 'horizontalalignment', 'center');     
        text(-rlen*sinct(i), -rlen * cosct(i), ticm2, 'horizontalalignment', 'center');
    end
    set(gca,'FontWeight', 'bold');

    axis('equal');
    axis('off');
    
existprn = unique(NMEAQM(:,5));
for aa = 1:length(existprn)
    prn = existprn(aa);
    if prn < 10
        PRN = strcat('G0',num2str(prn));
    elseif prn > 32
        PRN = strcat('GL',num2str(prn));
    else
        PRN = strcat('G',num2str(prn));
    end
    indxprn = find(NMEAQM(:,5) == prn);
    QM = NMEAQM(indxprn,:);
    Az = QM(:,7);     El = QM(:,6);
    yy = (El-90).* -(cos(Az*pi/180));
    xx = (El-90).* -(sin(Az*pi/180));
    ylast = yy(length(yy));
    xlast = xx(length(xx));
    figure(100);
    hold on;
%     plot(xx, yy, '.', 'Markersize', 15);
    text(xlast, ylast, PRN,'Fontsize',15)
end

% figure(100);
% hold on;
% yy = (El-90).* -(cos(Az*pi/180));
% xx = (El-90).* -(sin(Az*pi/180));
% 
% % 나타내고자하는 방위각, 고도각의 위치에 빨간색 점을 찍는다.
% plot(xx, yy, '.','color', 'r', 'Markersize', 15);
