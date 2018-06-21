function DDPlotQMcheck(QMfileBs, QMfileRv, type,START,STOP,event_time)
%
% function DDPlotQMcheck(QMfileBs, QMfileRv, type,START,STOP,event_time)
%
%   Read the Base, Rover QM and Plot 'QM'
%   
%   input QM1 : Base QM
%   input QM2 : Rover QM
%   input type : type of Observables
%               (111(L1), 112(L2), 120(C1), 121(P1), 122(P2), 123(C2), 
%               131(D1), 132(D2), 141(S1), 142(S2))
%
%   PRN = Figure 넘버링(입력 type 별)
%
%
%   Originally coded by Joonseong Gim, Mar 10, 2016
%
% clear all; close all;
% 
% QMfileBs = 'QSBBs16055';
% QMfileRv = 'QSBRv16056';
% type = 141;
%% Load QMfile(Base, Rover)
QM1 = load(QMfileBs);
QM2 = load(QMfileRv);

%% Base QMfile 내에 존재하는 PRN, Type
Base_prn_view = unique(QM1(:,2))';              % Base QMfile 내에 존재하는 PRN
disp('PRNs in Base QM')
disp(Base_prn_view)
Base_types = unique(QM1(:,3))';
disp('TYPEs in Base QM')
disp(Base_types)
%% Rover QMfile 내에 존재하는 PRN, Type
Rover_prn_view = unique(QM2(:,2))';             % Rover QMfile 내에 존재하는 PRN
disp('PRNs in Rover QM')
disp(Rover_prn_view)
Rover_types = unique(QM1(:,3))';
disp('TYPEs in Rover QM')
disp(Rover_types)

%% Type에 따른 Base, Rover QMfile Sorting
QM1 = SelectQM(QM1, type);
QM2 = SelectQM(QM2, type);

%% Type에 따른 Title 설정
switch type
    case 111
        TYPE = char('L1');
    case 112
        TYPE = char('L2');
    case 120
        TYPE = char('C1');
        MAX = 26000 ; MIN=20000;
        QM1(:,4) = QM1(:,4).*0.001;
        QM2(:,4) = QM2(:,4).*0.001
    case 121
        TYPE = char('P1');
    case 122
        TYPE = char('P2');
    case 123
        TYPE = char('C2');
    case 131
        TYPE = char('D1');
        MAX = 4000; MIN = -4000;
    case 132
        TYPE = char('D2');
    case 141
        TYPE = char('S1');
        MAX = 70; MIN=10;
    case 142
        TYPE = char('S2');
end


%% Base와 Rover 의 공통시간 추출
if START == 0 && STOP == 0
    FinalTTs = intersect(QM1(:, 1), QM2(:, 1));
else
    FinalTTs = [START:STOP];
end
%% 공통시간에 맞는 QM1, QM2 Sorting
QM1sort = QM1(min(find(QM1(:,1) == min(FinalTTs))):max(find(QM1(:,1) == max(FinalTTs))),:);
QM2sort = QM2(min(find(QM2(:,1) == min(FinalTTs))):max(find(QM2(:,1) == max(FinalTTs))),:);
%% 공통시간 내에 공통 PRN
QM1PRN = unique(QM1sort(:,2));
QM2PRN = unique(QM2sort(:,2));
commonprn = intersect(QM1PRN, QM2PRN);
EVENT = event_time(find(event_time(:,6) == 1),5);

for i = 1:length(commonprn)
    prn = commonprn(i);
%     prn = 28;
    PRN = num2str(prn);
    indxprn1 = find(QM1sort(:,2) == prn);
    indxprn2 = find(QM2sort(:,2) == prn);
    QM1prn = QM1sort(indxprn1,:);
    QM2prn = QM2sort(indxprn2,:);
    sametime = intersect(QM1prn(:,1), QM2prn(:,1));
    for j = 1:length(sametime)
        value1 = QM1prn(find(QM1prn(:,1) == sametime(j)),4);
        value2 = QM2prn(find(QM2prn(:,1) == sametime(j)),4);
        diff(j,1) = sametime(j);
        diff(j,2) = value1 - value2;
    end
    figure(prn)
    subplot(2,1,1)
    hold on; grid on;
    plot(QM1prn(:,1), QM1prn(:,4),'r.:')
    plot(QM2prn(:,1), QM2prn(:,4),'b.:')
    for ii = 1:length(EVENT)
        line([EVENT(ii) EVENT(ii)],[-9999 9999])
    end

    xlim([min(FinalTTs) max(FinalTTs)])
%     xlim([29570 29590]);
    ylim([MIN MAX])
    legend('Base','Rover')
    title(['PRN = ',PRN,'  TYPE = ',TYPE])
    xlabel('GPS second')
    
    subplot(2,1,2)
    plot(diff(:,1), diff(:,2),'r.:');
    for ii = 1:length(EVENT)
        line([EVENT(ii) EVENT(ii)],[-9999 9999])
    end

    xlim([min(FinalTTs) max(FinalTTs)])
%     xlim([29570 29590]);
    if type == 141;
        ylim([-25 25])
    elseif type == 131;
        ylim([-500 500])
    else
    end
    grid on;
%     ylabel(Y)
sametime = [];
diff = [];
end
