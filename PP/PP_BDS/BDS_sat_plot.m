clear all; close all

% navfile = 'brdm0250.17p';
% navfile = 'brdm0100.15p';

%% 연도를 선책하여 plot
%% yy=1(15년도), yy=2(17년도), yy=3(15,17동시)
yy = 3;
%% 시간간격 설정
tt = 300;
%% plot 저장 파일 포멧 설정
% filetype = '.fig';
filetype = '.png';

%% Beidou Leap second
LeapSecBDS = 14;

%% True Position
% TruePos = [-3041235.578 4053941.677 3859881.013];       % A point
% TruePos = [-3041241.741 4053944.143 3859873.640];       % B point
TruePos = [-2279828.8529 5004706.5404 3219777.4631]; %: JFNG site log 121029

%% Ephemeris 행렬 생성
% [eph, trashPRN, trashT]=ReadEPH_all(navfile);

%% Plot 과정 시작
%% 저장해 놓은 eph load
load('brdm010015p.mat');

%% Beidou Eph만 분리
eph = eph(find(eph(:,18) > 200 & eph(:,18) < 300),:);
prns1 = unique(eph(:,18));

%% 계산할 시간 설정
FinalTTs = [eph(1,8):tt:eph(1,8)+86400];
tHours = mod(FinalTTs, 86400);
tHours = tHours/3600;
cnt = 1;
for i=1:length(prns1)
    prn = prns1(i);
    for j=1:length(FinalTTs)-1
        gs = FinalTTs(j);
        tHour = tHours(j);
        icol = PickEPH(eph, prn, gs);
        toe = eph(icol, 8); a = eph(icol, 19); b = eph(icol, 20); c = eph(icol, 21); Tgd = eph(icol, 23);
        health = eph(icol, 22);
        STT = GetSTTbrdm(gs, eph, icol, TruePos(1:3)'); % 신호전달시간 계산
        tc = gs - STT- LeapSecBDS;  % 신호전달 시간 보정... bds윤초 반영
        SatPos = GetSatPosNC_GC(eph, icol, tc); % 위성위치 산출
        SatPos = RotSatPos(SatPos, STT); % 지구자전효과 고려
        fprintf('%3d : %6d\n', prn, gs);
        results1(cnt,:) = [prn, gs, tHour, SatPos', health];
        cnt = cnt + 1;
    end
end

load('brdm025017p.mat');

%% Beidou Eph만 분리
eph = eph(find(eph(:,18) > 200 & eph(:,18) < 300),:);
prns2 = unique(eph(:,18));

%% 계산할 시간 설정
FinalTTs = [eph(1,8):tt:eph(1,8)+86400];
tHours = mod(FinalTTs, 86400);
tHours = tHours/3600;
cnt = 1;
for i=1:length(prns2)
    prn = prns2(i);
    for j=1:length(FinalTTs)-1
        gs = FinalTTs(j);
        tHour = tHours(j);
        icol = PickEPH(eph, prn, gs);
        toe = eph(icol, 8); a = eph(icol, 19); b = eph(icol, 20); c = eph(icol, 21); Tgd = eph(icol, 23);
        health = eph(icol, 22);
        STT = GetSTTbrdm(gs, eph, icol, TruePos(1:3)'); % 신호전달시간 계산
        tc = gs - STT- LeapSecBDS;  % 신호전달 시간 보정... bds윤초 반영
        SatPos = GetSatPosNC_GC(eph, icol, tc); % 위성위치 산출
        SatPos = RotSatPos(SatPos, STT); % 지구자전효과 고려
        fprintf('%3d : %6d\n', prn, gs);
        results2(cnt,:) = [prn, gs, tHour, SatPos', health];
        cnt = cnt + 1;
    end
end

totalprn = unique([prns1;prns2]);
switch yy
    case 1
        for i=1:length(prns1)
            prn = prns1(i);
            PRN = num2str(prns1(i));
            figure(prn)
            subplot(4,1,1)
            suptitle(['PRN : ',PRN,'(15)'])
            hold on; grid on;
            plot(results1(find(results1(:,1) == prn),3),results1(find(results1(:,1) == prn),4),'r.:')
            ylabel('Sat Pos x')
            xlim([0,24])
            subplot(4,1,2)
            hold on; grid on;
            plot(results1(find(results1(:,1) == prn),3),results1(find(results1(:,1) == prn),5),'r.:')
            ylabel('Sat Pos y')
            xlim([0,24])
            subplot(4,1,3)
            hold on; grid on;
            plot(results1(find(results1(:,1) == prn),3),results1(find(results1(:,1) == prn),6),'r.:')
            ylabel('Sat Pos z')
            xlim([0,24])
            subplot(4,1,4)
            hold on; grid on;
            stairs(results1(find(results1(:,1) == prn),3),results1(find(results1(:,1) == prn),7),'ko-')
            xlabel('hours')
            ylabel('Sat Health')
            xlim([0,24])
            figname = strcat('15prn',PRN,'_',num2str(tt),filetype);
%             saveas(gcf,figname)
        end
    case 2
        for i=1:length(prns2)
            prn = prns2(i);
            PRN = num2str(prns2(i));
            figure(prn)
            subplot(4,1,1)
            suptitle(['PRN : ',PRN,'(17)'])
            hold on; grid on;
            plot(results2(find(results2(:,1) == prn),3),results2(find(results2(:,1) == prn),4),'b.:')
            ylabel('Sat Pos x')
            xlim([0,24])
            subplot(4,1,2)
            hold on; grid on;
            plot(results2(find(results2(:,1) == prn),3),results2(find(results2(:,1) == prn),5),'b.:')
            ylabel('Sat Pos y')
            xlim([0,24])
            subplot(4,1,3)
            hold on; grid on;
            plot(results2(find(results2(:,1) == prn),3),results2(find(results2(:,1) == prn),6),'b.:')
            ylabel('Sat Pos z')
            xlim([0,24])
            subplot(4,1,4)
            hold on; grid on;
            stairs(results2(find(results2(:,1) == prn),3),results2(find(results2(:,1) == prn),7),'ko-')
            xlabel('hours')
            ylabel('Sat Health')
            xlim([0,24])
            figname = strcat('17prn',PRN,'_',num2str(tt),filetype);
%             saveas(gcf,figname)
        end
    case 3
        for i=1:length(totalprn)
            prn = totalprn(i);
            PRN = num2str(totalprn(i));
            figure(prn)
            subplot(4,1,1)
            suptitle(['PRN : ',PRN,'(15(red) & 17(blue)) '])
            hold on; grid on;
            if ~isempty(find(prns1(:,1) == prn)) & ~isempty(find(prns2(:,1) == prn))
                plot(results1(find(results1(:,1) == prn),3),results1(find(results1(:,1) == prn),4),'r.:')
                plot(results2(find(results2(:,1) == prn),3),results2(find(results2(:,1) == prn),4),'b.:')
                
            elseif ~isempty(find(prns1(:,1) == prn))
                plot(results1(find(results1(:,1) == prn),3),results1(find(results1(:,1) == prn),4),'r.:')
            elseif ~isempty(find(prns2(:,1) == prn))
                plot(results2(find(results2(:,1) == prn),3),results2(find(results2(:,1) == prn),4),'b.:')
            end
            ylabel('Sat Pos x')
            xlim([0,24])
            
            subplot(4,1,2)
            hold on; grid on;
            if ~isempty(find(prns1(:,1) == prn)) & ~isempty(find(prns2(:,1) == prn))
                plot(results1(find(results1(:,1) == prn),3),results1(find(results1(:,1) == prn),5),'r.:')
                plot(results2(find(results2(:,1) == prn),3),results2(find(results2(:,1) == prn),5),'b.:')
                
            elseif ~isempty(find(prns1(:,1) == prn))
                plot(results1(find(results1(:,1) == prn),3),results1(find(results1(:,1) == prn),5),'r.:')
            elseif ~isempty(find(prns2(:,1) == prn))
                plot(results2(find(results2(:,1) == prn),3),results2(find(results2(:,1) == prn),5),'b.:')
            end
            ylabel('Sat Pos y')
            xlim([0,24])
            
            subplot(4,1,3)
            hold on; grid on;
            if ~isempty(find(prns1(:,1) == prn)) & ~isempty(find(prns2(:,1) == prn))
                plot(results1(find(results1(:,1) == prn),3),results1(find(results1(:,1) == prn),6),'r.:')
                plot(results2(find(results2(:,1) == prn),3),results2(find(results2(:,1) == prn),6),'b.:')
                
            elseif ~isempty(find(prns1(:,1) == prn))
                plot(results1(find(results1(:,1) == prn),3),results1(find(results1(:,1) == prn),6),'r.:')
            elseif ~isempty(find(prns2(:,1) == prn))
                plot(results2(find(results2(:,1) == prn),3),results2(find(results2(:,1) == prn),6),'b.:')
            end
            ylabel('Sat Pos z')
            xlim([0,24])
            
            subplot(4,1,4)
            hold on; grid on;
            if ~isempty(find(prns1(:,1) == prn)) & ~isempty(find(prns2(:,1) == prn))
                stairs(results1(find(results1(:,1) == prn),3),results1(find(results1(:,1) == prn),7),'r.:')
                stairs(results2(find(results2(:,1) == prn),3),results2(find(results2(:,1) == prn),7),'b.:')
                
            elseif ~isempty(find(prns1(:,1) == prn))
                stairs(results1(find(results1(:,1) == prn),3),results1(find(results1(:,1) == prn),7),'r.:')
            elseif ~isempty(find(prns2(:,1) == prn))
                stairs(results2(find(results2(:,1) == prn),3),results2(find(results2(:,1) == prn),7),'b.:')
            end
            xlabel('hours')
            ylabel('Sat Health')
            xlim([0,24])
            figname = strcat('1517prn',PRN,'_',num2str(tt),filetype);
%             saveas(gcf,figname)
        end
end

