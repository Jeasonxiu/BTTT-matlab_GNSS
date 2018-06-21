function [QM_sorted] = PlotQM(QMfile, prn, type)
%
%function [QM] = PlotQM(QMfile, prn, type)
%
%   Read the given QMfile and Plot 'QM_sorted'
%   
%   input QMfile : WriteObs
%   input prn : prn number
%   input type : type of Observables
%               (111(L1), 112(L2), 120(C1), 121(P1), 122(P2), 123(C2), 
%               131(D1), 132(D2), 141(S1), 142(S2))
%
%   Output QM : QM Matrix[gs prn type obs]
%
%   Example : result = get_eph('test.txt', svn, type)
%
%   Originally coded by Joonseong Gim, Jan 5, 2016
%

% prn = 2;
% type = 120;
% load('QMfile');

%% QMfile 내 존재하는 prn, type 정의
prn_view = unique(QMfile(:,2))';
disp('PRNs in QMfile')
disp(prn_view)
types = unique(QMfile(:,3))';
disp('TYPEs in QMfile')
disp(types)
gs = unique(QMfile(:,1));

%% 선택된 prn, type에 따른 QM_Sorted matrix 생성
ObsType = type;
if isempty(find(prn_view==prn))
    sprintf('prn not exist')
    QM_sorted = 0;
    
elseif isempty(find(types==type))
    sprintf('prn exist, but type not exist')
    
else
    [arrQM, FinalPRNs, FinalTTs] = ReadQM('QMfile'); 
    QM = SelectQM(arrQM, ObsType);
    indxprn = find(QM(:,2) == prn);
    if length(indxprn) == 0
        return
    else
        QM = QM(indxprn,:);
    end
    
    switch type
        case 111
            TYPE = char('L1');
        case 112
            TYPE = char('L2');
        case 120
            TYPE = char('C1');
        case 121
            TYPE = char('P1');
        case 122
            TYPE = char('P2');
        case 123
            TYPE = char('C2');
        case 131
            TYPE = char('D1');
        case 132
            TYPE = char('D2');
        case 141
            TYPE = char('S1');
        case 142
            TYPE = char('S2');
    end
            
figure(101)            
plot(QM(:,1),QM(:,4),'r.:')
grid on
xlim([min(QM(:,1)) max(QM(:,1))])
title({'\fontsize{16}PlotQM';['\fontsize{11}PRN ',num2str(prn),', ','TYPE ',num2str(type),'(',TYPE,')']})
xlabel('\fontsize{12}gs')
ylabel('\fontsize{12}obs value')    
 
QM_sorted = QM;
end

for i = 1:length(commonprn)
    prn = commonprn(i);
%     prn = 28;
    PRN = num2str(prn);
    indxprn1 = find(QM1sort(:,2) == prn);
    indxprn2 = find(QM2sort(:,2) == prn);
    QM1prn = QM1sort(indxprn1,:);
    QM2prn = QM2sort(indxprn2,:);
    sametime = intersect(QM1prn(:,1), QM2prn(:,1));
    if ~isempty(sametime)
        for j = 1:length(sametime)
            value1 = QM1prn(find(QM1prn(:,1) == sametime(j)),4);
            value2 = QM2prn(find(QM2prn(:,1) == sametime(j)),4);
            diff(j,1) = sametime(j);
            diff(j,2) = value1 - value2;
        end
    else
        diff = [];
    end
    figure(prn)
    subplot(2,1,1)
    hold on; grid on;
    plot(QM1prn(:,1), QM1prn(:,4),'r.:')
    plot(QM2prn(:,1), QM2prn(:,4),'b.:')
%     line([307782 307782],[-9999 9999])
%     line([307786 307786],[-9999 9999])
    xlim([min(FinalTTs) max(FinalTTs)])
%     xlim([29570 29590]);
    ylim([MIN MAX])
    legend(role1, role2)
    title(['PRN = ',PRN,'  TYPE = ',TYPE])
    xlabel('GPS second')
    subplot(2,1,2)
    if ~isempty(diff)
        plot(diff(:,1), diff(:,2),'r.:');
        legend(strcat(role1, ' SNR - ', role2, ' SNR'))
%         line([307782 307782],[-9999 9999])
%         line([307786 307786],[-9999 9999])
        xlim([min(FinalTTs) max(FinalTTs)])
        %     xlim([29570 29590]);
        ylim([-20 20])
        grid on;
    else
    end
        %     ylabel(Y)
sametime = [];
diff = [];
end







