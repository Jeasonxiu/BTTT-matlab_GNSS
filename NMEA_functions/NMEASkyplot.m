% function NMEASkyplot(NMEAfile)
%=====================================================
% DO: After reading QM file and then draw skyplot
% Input  : NMEAfile
% Output : figure or jpg file

% clear all; close all;
% NMEAfile = 'jprA055f_vu.txt';



% POINT = 'A';                    % logginf Point
% DOY = '055'; YY = '16';         % Day of Year
% TIME = 'f';                     % Logging Hour(9a,10b,11c,12d,13e,14f,15g,16h,17i,18j,19k,
% %             20l,21m,22n,23o,24p,1q,2r,3s,4t,5u,6v,7w,8x)
% DEVICE = 'vu';
% NMEAfile = strcat('jpr',POINT,DOY,TIME,'_',DEVICE,'.txt');
% NMEAQM = writeNMEA(NMEAfile);

cutoff = 15;
site = 'E';

figure(100);
Skymap(cutoff);
% Az = NMEAQM(:,7);
% El = NMEAQM(:,6);
% figure(100);
% hold on;
% yy = (El-90).* -(cos(Az*pi/180));
% xx = (El-90).* -(sin(Az*pi/180));

% 나타내고자하는 방위각, 고도각의 위치에 빨간색 점을 찍는다.
% plot(xx, yy, '.','color', 'r', 'Markersize', 15);

existprn = unique(NMEAQM(:,5));
existprn = existprn(existprn(:,1) > 0, :);
for aa = 1:length(existprn)
    prn = existprn(aa);
    
    if prn < 10
        PRN = strcat('G0',num2str(prn));
    elseif prn >70
        PRN = strcat('GL',num2str(prn));
    else
        PRN = strcat('G',num2str(prn));
    end
    indxprn = find(NMEAQM(:,5) == prn);
    QM = NMEAQM(indxprn,:);
    indexzero = find(QM(:,6) ~= 0);
    QM = QM(indexzero,:); 
    indexzero = find(QM(:,7) ~= 0); 
    QM = QM(indexzero,:);    
    Az = QM(:,7);     El = QM(:,6);
    %     yy = (El-90).* -(cos(Az*pi/180));
    %     xx = (El-90).* -(sin(Az*pi/180));
    yy = (El-90).* -(cos(Az*pi/180));
    xx = (El-90).* -(sin(Az*pi/180));
    ylast = yy(length(yy));
    xlast = xx(length(xx));
    figure(100);
    hold on;
    plot(xx, yy, '.', 'Markersize', 25);
    plot(xlast, ylast, 'r.', 'Markersize', 25);
    text(xlast, ylast, PRN,'Fontsize',15)
%     title('LG vu, site = A, DOY=055')
end