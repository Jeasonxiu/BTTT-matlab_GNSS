function PPPlotQM(QM1,type)
% 
% function PPPlotQM(QMfileBs, QMfileRv, type)
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

%% Load QMfile(Base, Rover)


%% QMfile 내에 존재하는 PRN, Type
prn_view = unique(QM1(:,2))';              % Base QMfile 내에 존재하는 PRN
disp('PRNs in QMfile')
disp(prn_view)
types = unique(QM1(:,3))';
disp('TYPEs in QMfile')
disp(types)

%% Type에 따른 Title 설정
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
        MAX = 4000; MIN = -4000;
    case 132
        TYPE = char('D2');
    case 141
        TYPE = char('S1');
        MAX = 70; MIN=10;
    case 142
        TYPE = char('S2');
end

%% Type에 따른 Base, Rover QMfile Sorting
% QM1 = SelectQM(QM1, type);

commonprn = unique(QM1(:,2));

for i = 1:length(commonprn)
    prn = commonprn(i);
%     prn = 28;
    PRN = num2str(prn);
    indxprn1 = find(QM1(:,2) == prn);
    QM1prn = QM1(indxprn1,:);
    
    figure(prn)
    hold on; grid on;
    plot(QM1prn(:,1), QM1prn(:,4),'r.:')
%     line([307782 307782],[-9999 9999])
%     line([307786 307786],[-9999 9999])
    xlim([min(QM1(:,1)) max(QM1(:,1))])
    ylim([MIN MAX])
    title(['PRN = ',PRN,'  TYPE = ',TYPE])
    xlabel('GPS second')
    
end