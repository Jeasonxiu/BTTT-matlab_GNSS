% function SatPosDiff(SatPos)
clear all
close all

load('cp_alm_test.mat');
finalprn = unique(SatPos(:,2));
finalTTs = unique(SatPos(:,1));
tHour = mod(finalTTs, 86400);
tHour = tHour/3600;

for i = 1: length(finalprn)
    prn = finalprn(i);
    Sat = []; Diff = []; Diff_norm = [];
    Sat = SatPos(find(SatPos(:,2) == prn), :);
    Diff = [mod(Sat(:,1), 86400)/3600, Sat(:,3:5) - Sat(:,6:8)];
    for j = 1:length(Diff(:,1))
        Diff_norm(j,1:2) = [Diff(j,1), norm(Diff(j,2:4))];
    end
    if i < 6
        figure(1)
        subplot(5,1,[i])
        plot(Diff(:,1),Diff_norm(:,2),'b.');
        hold on; grid on;
        xlim([min(tHour) max(tHour)]);
        PRN = strcat('PRN ',num2str(prn));
        title(PRN)  
        ylabel('(meter)')
        xlabel('Hours');
    elseif i < 11
        figure(2)
        subplot(5,1,[i-5])
        plot(Diff(:,1),Diff_norm(:,2),'b.');
        hold on; grid on;
        xlim([min(tHour) max(tHour)]);
        PRN = strcat('PRN ',num2str(prn));
        title(PRN)  
        ylabel('(meter)')
        xlabel('Hours');
    elseif i < 16
        figure(3)
        subplot(5,1,[i-10])
        plot(Diff(:,1),Diff_norm(:,2),'b.');
        hold on; grid on;
        xlim([min(tHour) max(tHour)]);
        PRN = strcat('PRN ',num2str(prn));
        title(PRN)  
        ylabel('(meter)')
        xlabel('Hours');
    end
end
