clear all
close all

% QM = load('QDSAJ_18135');    % 18135 JAVAD
% QM = load('QDSBJ_18135');    % 18135 S8
% QM = load('QDSBJ_18137');    % 18137 S8
% QM = load('QDSBS_18145');    % 18145 S8
% QM = load('QDSBN_18145');    % 18145 Nexus9
% QM = load('QDSBS_18147');    % 18147 S8
% QM = load('QDSBS1_18148');    % 18148 S8 1
% QM = load('QUSN_16279');    % 16300 N9
% QM = load('QDSBS2_18148');    % 18148 S8 2
% QM = load('QDSBN2_18148');    % 18148 N9 2
% QM = load('QIHU4_18149');    % 18149 INHA Note8
% QM = load('QDSBNU_18149');    % 18149 NEXUS9 up
% QM = load('QDSBND_18149');    % 18149 NEXUS9 down
% QM = load('QSHNU_18150');    % 18150 NEXUS9 up 신항
% QM = load('QSHND_18150');    % 18150 NEXUS9 down 신항
QM = load('QSHS1_18150');    % 18150 S8 down 신항
% load('18165_rt_data.mat');
% QM = QDBN8_18165;    % 18165 S9(GEO++) L1 : cycle
QM(:,1) = round(QM(:,1));
QM_C1 = QM(QM(:,3) == 120, :);        % Base C1 데이터 추출
QM_L1 = QM(QM(:,3) == 111, :);        % Base L1 데이터 추출
% QM_L1(:,4) = QM_L1(:,4)/0.1903;

prnlist = unique(QM_C1(:,2));
FinalTTs = unique(QM(:,1));
% FinalTT = FinalTTs(1:515);          % 18150 NEXUS9 U 신항
% FinalTT = FinalTTs(516:963);        % 18150 NEXUS9 U 신항
% FinalTT = FinalTTs(964:1351);       % 18150 NEXUS9 U 신항
% FinalTT = FinalTTs(1352:1810);      % 18150 NEXUS9 U 신항
% FinalTT = FinalTTs(1811:2185);      % 18150 NEXUS9 U 신항
% FinalTT = FinalTTs(2186:end);       % 18150 NEXUS9 U 신항
% FinalTTs = FinalTT;      % 18150 NEXUS9 U 신항
tHour = mod(FinalTTs, 86400); tHour = tHour/3600;
if find(tHour(:) == 0) > 1
    tHour(find(tHour(:) == 0):end) = tHour(find(tHour(:) == 0):end) + 24;
end
        
Culoum = zeros(length(prnlist),200);
Culoum_gs = zeros(length(prnlist),200);
Culoum_diff = zeros(length(prnlist),200);
for i=1:length(prnlist)
    prn = prnlist(i);
    qm_l1 = QM_L1(find(QM_L1(:,2) == prn),:);
    prn_tHour = mod(qm_l1(:,1), 86400); prn_tHour = prn_tHour/3600;
    if find(prn_tHour(:) == 0) > 1
        prn_tHour(find(prn_tHour(:) == 0):end) = prn_tHour(find(prn_tHour(:) == 0):end) + 24;
    end
    l1_diff = diff(qm_l1);
    if ~isempty(find(abs(l1_diff(:,4)) > 1000)')
        Culoum(i,1:length(find(abs(l1_diff(:,4)) > 1000))+1) =...
            [prn, find(abs(l1_diff(:,4)) > 1000)'+1];
        Culoum_diff(i,1:length(find(abs(l1_diff(:,4)) > 1000))+1) =...
            [prn, l1_diff(find(abs(l1_diff(:,4)) > 1000)',4)'];
        Culoum_gs(i,1:length(find(abs(l1_diff(:,4)) > 1000))+1) =...
            [prn, qm_l1(find(abs(l1_diff(:,4)) > 1000)'+1,1)'];
    end
    figure(prn)
    subplot(3,1,1)
    plot(prn_tHour, qm_l1(:,4),'b.')
    grid on; hold on;
    xlim([min(tHour), max(tHour)]);
    ylabel('Meters')
    subplot(3,1,2)
    plot(prn_tHour(2:end), diff(qm_l1(:,4)),'b.-')
    grid on; hold on;
    xlim([min(tHour), max(tHour)]);
    subplot(3,1,3)
    plot(prn_tHour(3:end), diff(qm_l1(:,4),2),'b.-')
    grid on; hold on;
    xlim([min(tHour), max(tHour)]);
    ylim([-3, 3])
    xlabel('Hour')
end


tHour = mod(QM_L1(:,1), 86400); tHour = tHour/3600;
if find(tHour(:) == 0) > 1
    tHour(find(tHour(:) == 0):end) = tHour(find(tHour(:) == 0):end) + 24;
end
QM_L1(:,1) = tHour;
for i=1:length(tHour)
    gs = tHour(i);
    Epoch = QM_L1(find(QM_L1(:,1) == gs),:);
    
    numprn(i,:) = [gs, length(Epoch(:,2))];
end
figure(99)
stairs(numprn(:,1), numprn(:,2),'bo-')
grid on; 
xlim([min(numprn(:,1)), max(numprn(:,1))]);
xlabel('Hour')
ylim([1, 12])
    
mean(numprn(:,2))
% figure(199)
% stairs(states(:,1), states(:,2),'bo-')
% grid on; 
% xlim([min(states(:,1)), max(states(:,1))]);
% xlabel('Hour')
% ylim([0, 5])