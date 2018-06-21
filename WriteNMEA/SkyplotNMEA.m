function SkyplotNMEA(NMEAQM)
% SYNTAX:
%   SkyplotNMEA(NMEAQM)
%
% INPUT:
%   NMEAQM : after writeNMEA(filename, year, month, day)
%
% OUTPUT:
%   SkyPlot

cutoff = 15;

figure(100);
Skymap(cutoff);

existprn = unique(NMEAQM(:,7));
existprn = existprn(existprn(:,1) > 0, :);
for aa = 1:length(existprn)
    prn = existprn(aa);
    if prn < 110
        PRN = strcat('G0',num2str(prn));
        flag = 1;
    elseif prn > 109 & prn < 200
        PRN = strcat('G',num2str(prn));
        flag = 1;
    elseif prn > 300 & prn < 400
        PRN = strcat('R',num2str(prn));
        flag = 2;
    elseif prn > 200 & prn < 300
        PRN = strcat('C',num2str(prn));
        flag = 3;
    elseif prn > 500
        PRN = strcat('QZ',num2str(prn));
        flag = 4;
    end
    indxprn = find(NMEAQM(:,7) == prn);
    QM = NMEAQM(indxprn,:);
    indexzero = find(QM(:,8) ~= 0);
    QM = QM(indexzero,:);
    indexzero = find(QM(:,9) ~= 0);
    QM = QM(indexzero,:);
    if ~isempty(QM)
        Az = QM(:,9);     El = QM(:,8);
        yy = (El-90).* -(cos(Az*pi/180));
        xx = (El-90).* -(sin(Az*pi/180));
        ylast = yy(length(yy));
        xlast = xx(length(xx));
        figure(100);
        hold on;
        switch flag
            case 1
                plot(xx, yy, '.', 'Markersize', 15);
                plot(xlast, ylast, 'r.', 'Markersize', 15);
            case 2
                plot(xx, yy, '*', 'Markersize', 5);
                plot(xlast, ylast, 'r*', 'Markersize', 5);
            case 3
                plot(xx, yy, 'd', 'Markersize', 5);
                plot(xlast, ylast, 'rd', 'Markersize', 5);
            case 1
                plot(xx, yy, '+', 'Markersize', 5);
                plot(xlast, ylast, 'r+', 'Markersize', 5);
        end
        text(xlast, ylast, PRN,'Fontsize',15)
    end
    flag = 0;
end