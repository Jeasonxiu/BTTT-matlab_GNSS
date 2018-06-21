% function SkyPlotQM(QM, x, year, DOY, cutoff, interval, system)
%
%   SkyPlotQM(QMfile,DOY,year,interval)
%   : QM array에 존재하는 위성항법시스템의 해당시간 Skyplot.
%     각 위성항법 시스템별 Skyplot(N 개), 통한 Skyplot
%
%   input 1) QMfile   = QM array
%         2) x        = vec_site(1X3) : ecef
%         3) year     = year : for navfile [YY]
%         4) DOY      = day of year : for navfile [DOY]
%         5) cutoff   = Cutoff elevation angle [XX]
%         6) interval = Satellite Pos interval [ss]
%         7) system   = select system
%                       a : all, g : GPS, r : GLONASS, c : Beidou
%
% Copyright: Joonseong Gim, 6/03/2018
%
close all;
clear all; 
% %% test parameters
% load('rover_B_1.mat');
% QM = load('QTHu1_18059');        % 2018-02-28 teheran javad1
% QM = load('QDSBN_18151');        % 2018-05-31 Note8
% QM = load('QDSBN_18152');        % 2018-06-01 Note8
% QM = load('QDSBN_18156');        % 2018-06-05 Note8
% QM = load('QDSBS_18156');        % 2018-06-05 Note8
% QM = load('QDSBN4_18156');        % 2018-06-05 Note8
% load('eph180605.mat');
% QM = load('QIHUN_18157');       % 2018-06-06 Note8 Inha 30분
% QM = load('QIHUN_18157_whole');       % 2018-06-06 Note8 Inha 
QM = load('QDSBS1_18157');       % 2018-06-06 Note8 Inha 
load('eph180606.mat');
QM(:,1) = round(QM(:,1));
x = [-3053196.72849527,4042608.64018026,3862339.72264547];
% DOY = 059; year = 18;
% DOY = 151; year = 18;           % 180531
% DOY = 152; year = 18;           % 180601
% DOY = 156; year = 18;           % 180605
DOY = 157; year = 18;           % 180605
cutoff = 15;
interval = 30;
system ='g';
% load('eph180228.mat');      % 2018-02-28 teheran

%% start
x_gd = xyz2gd(x);       % 현재 위치
[gw, gd] = ydoy2gwgd(year, DOY); %: GPS WEEK 결정

%% 윤초 핸들링
LeapSecBDS = 14;        % BDS
LeapSecGLO = 18;           % GLO

%% 측정시간 체크 및 QM 핸들링
FinalTTs = unique(QM(:,1));
% QM(find(QM(:,3) < 200),2) = QM(find(QM(:,3) < 200),2) + 100;      % GPS 위성 +100
% QM(find(QM(:,3) > 300),2) = QM(find(QM(:,3) > 300),2) + 300;      % GLO 위성 +300
% QM(find(QM(:,3) < 300 & QM(:,3) > 200),2) = QM(find(QM(:,3) < 300 & QM(:,3) > 200),2) + 00;      % BDS 위성 +200
SNRinQM = QM(find(QM(:,3) == 141 | QM(:,3) == 241 | QM(:,3) == 341), :);    % SNR 관측값 추출
SNRinQM(:,4) = round(SNRinQM(:,4));
%% Navigation file load
FileNav = strcat('brdm', num2str(DOY,'%03d'), '0.', num2str(year,'%02d'), 'p');   %: Navigation RINEX file
disp(['reading a ', FileNav, ' nav file!'])
tic;
% [eph, trashPRN, trashT]=ReadEPH_all(FileNav);     % make a eph array
toc;
disp('eph array making complete.')

%% GLONASS eph 생성
ephGLO = eph;
ephGLO=ephGLO(ephGLO(:,1)<400&ephGLO(:,1)>300,:);ephGLO(:,1)=ephGLO(:,1);

%% SV pos : 현재 관측된 위성, 항법메시지 위성
PRNineph = unique(eph(find(eph(:,18) < 400 & eph(:,18) ~= 104), 18));      % eph 어레이에 있는 PRN
%% system 선택
system = upper(system);
switch system
    case 'G'
%         QM = QM(find(QM(:,3) < 200),:);
        PRNineph = PRNineph(find(PRNineph < 200));
    case 'R'
%         QM = QM(find(QM(:,3) > 300 & QM(:,4) < 400),:);
        PRNineph = PRNineph(find(PRNineph > 300 & PRNineph < 400));
    case 'C'
%         QM = QM(find(QM(:,3) > 200 & QM(:,4) < 300),:);
        PRNineph = PRNineph(find(PRNineph > 200 & PRNineph < 300));
end
datafromeph = [];
disp('Sat Pos processing!')
tic;
for i = 1:interval:length(FinalTTs)
    gs = FinalTTs(i);
    PRNinQM1e = QM(find(QM(:,1) == gs),2);
    for j = 1:length(PRNineph)
        prn = PRNineph(j);
        if ~isempty(find(PRNinQM1e == prn))
            snr = SNRinQM(find(SNRinQM(:,1) == gs & SNRinQM(:,2) == prn),4);
        else
            snr = 0;
        end
        
        if prn < 300
            icol = PickEPH(eph, prn, gs);
            STT = GetSTTbrdm(gs, eph, icol, x(1:3)'); % 신호전달시간 계산
            if prn < 200
                tc = gs - STT ;             % 신호전달시간 보정
            else
                tc = gs - STT- LeapSecBDS;  % 신호전달 시간 보정... bds윤초 반영
            end
            SatPos = GetSatPosNC_GC(eph, icol, tc); % 위성위치 산출
            SatPos = RotSatPos(SatPos, STT); % 지구자전효과 고려
            
        elseif prn > 300
            icol = PickEPH_GLO2(ephGLO, prn, gs);
            STT = GetSTTbrdc_GLO(gs, ephGLO, icol, x(1:3)', 60); % 방송궤도력으로 신호전달시간 계산
            tc = gs - STT;
            tc = tc - LeapSecGLO;
            [SatPos,SatVel,SatLS] = GetSatPosGLO(ephGLO, icol, tc, 60); % 방송궤도력으로 위성위치 계산
            SatPos = RotSatPos(SatPos, STT); % 지구자전효과 고려
        end
        vec_rho = SatPos - x(1:3)';
        [az,el] = xyz2azel(vec_rho, x_gd(1), x_gd(2));
        SV(j,:) = [gs, prn, SatPos', az, el, snr];
    end
    datafromeph = [datafromeph; SV];      % Sat's pos, azimuth, elevation
end
toc;
disp('Sat Pos processing complete.')
datafromeph = datafromeph(find(datafromeph(:,7) > cutoff),:);
datafromQM = datafromeph(find(datafromeph(:,8) ~= 0),:);

%% eph skyplot
Az = datafromeph(:,6); El = datafromeph(:,7);
yy = (El-90).* -(cos(Az*pi/180));
xx = (El-90).* -(sin(Az*pi/180));
datafromeph(:,9) = xx; datafromeph(:,10) = yy;

figure(1)
hold on; grid on;
Skymap(cutoff);
[yys, mos, dds, hhs, mms, sss] = gwgs2date(gw, min(datafromeph(:,1)));
[yye, moe, dde, hhe, mme, sse] = gwgs2date(gw, max(datafromeph(:,1)));
title([num2str(yys),'/',num2str(mos),'/',num2str(dds),'  ',num2str(hhs),':',num2str(mms),':',num2str(sss),' ~ ',...
    num2str(yye),'/',num2str(moe),'/',num2str(dde),'  ',num2str(hhe),':',num2str(mme),':',num2str(sse),'(UTC)'])

prn = unique(datafromeph(:,2))';
disp('** Figure 1 is brdc SkyPlot**')
for i=1:length(prn)
    prn_ = prn(i);
    
    data_ = datafromeph(datafromeph(:,2) == prn_, :);
    len = length(data_(:,1));
    if prn_ < 110
        PRN = strcat('G0',num2str(prn_-100));
    elseif prn_ < 200
        PRN = strcat('G ',num2str(prn_-100));
    elseif prn_ > 200 & prn_ < 210
        PRN = strcat('C0',num2str(prn_-200));
    elseif prn_ > 209 & prn_ < 300
        PRN = strcat('C ',num2str(prn_-200));
    elseif prn_ > 300 & prn_ < 310
        PRN = strcat('R0',num2str(prn_-300));
    elseif prn_ > 309 & prn_ < 400
        PRN = strcat('R ',num2str(prn_-300));
    end
    
    xlast = data_(len,9);       % 위성 마지막 위치
    ylast = data_(len,10);      % 위성 마지막 위치
    if prn_ < 200
        for m = 1:len
            x_ = data_(m,9); y_ = data_(m,10);
            plot(x_, y_, 'r.','Markersize',10);
            text(xlast, ylast, PRN,'Fontsize',10)
        end
    elseif prn_ > 200 & prn_ < 300
        for m = 1:len
            x_ = data_(m,9); y_ = data_(m,10);
            plot(x_, y_, 'b.','Markersize',10);
            text(xlast, ylast, PRN,'Fontsize',10)
        end
    elseif prn_ > 300
        for m = 1:len
            x_ = data_(m,9); y_ = data_(m,10);
            plot(x_, y_, 'g.','Markersize',10);
            text(xlast, ylast, PRN,'Fontsize',10)
        end
    end
end

%% QM skyplot
Az = datafromQM(:,6); El = datafromQM(:,7);
yy = (El-90).* -(cos(Az*pi/180));
xx = (El-90).* -(sin(Az*pi/180));
datafromQM(:,9) = xx; datafromQM(:,10) = yy;

figure(2)
hold on; grid on;
Skymap(cutoff);
title([num2str(yys),'/',num2str(mos),'/',num2str(dds),'  ',num2str(hhs),':',num2str(mms),':',num2str(sss),' ~ ',...
    num2str(yye),'/',num2str(moe),'/',num2str(dde),'  ',num2str(hhe),':',num2str(mme),':',num2str(sse),'(UTC)'])
if length(datafromQM(1,:)) > 7
    map = jet(60);
    colormap(map)
end

prn = unique(datafromQM(:,2))';
disp('** Figure 2 is QM SkyPlot**')
for i=1:length(prn)
    prn_ = prn(i);
    
    data_ = datafromQM(datafromQM(:,2) == prn_, :);
    len = length(data_(:,1));
    if prn_ < 110
        PRN = strcat('G0',num2str(prn_-100));
    elseif prn_ < 200
        PRN = strcat('G ',num2str(prn_-100));
    elseif prn_ > 200 & prn_ < 210
        PRN = strcat('C0',num2str(prn_-200));
    elseif prn_ > 209 & prn_ < 300
        PRN = strcat('C ',num2str(prn_-200));
    elseif prn_ > 300 & prn_ < 310
        PRN = strcat('R0',num2str(prn_-300));
    elseif prn_ > 309 & prn_ < 400
        PRN = strcat('R ',num2str(prn_-300));
    end
    
    xlast = data_(len,9);       % 위성 마지막 위치
    ylast = data_(len,10);      % 위성 마지막 위치
    
    for m = 1:len
        color = map(data_(m,8),:);
        x_ = data_(m,9); y_ = data_(m,10);
        plot(x_, y_, '.','Markersize',10,'color', color);
%         plot(x_, y_, '*','color', color);
        text(xlast, ylast, PRN,'Fontsize',10)
    end
end

if length(datafromQM(1,:)) > 7
    colorbar('Ticks',[0.25, 0.5, 0.75],...
        'Ticklabels',{'30', '40', '50'});
end
%
%
%
%
%
