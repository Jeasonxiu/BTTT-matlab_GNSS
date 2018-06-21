clear all
close all

% load('18165_rt_data.mat');
% QM = QDBN8_18165;    % 18165 S9(GEO++) L1 : cycle
% QM = load('QDSBN_18151');    % 18151 Note8 Whole
% QM = load('QDSBN1_18151');    % 18151 Note8 B-1
% QM = load('QDSBN2_18151');    % 18151 Note8 B-2
% QM = load('QDSBN3_18151');    % 18151 Note8 B-3
% QM = load('QDSBN4_18151');    % 18151 Note8 B-4
% QM = load('QDSBN_18152');    % 18152 Note8 Whole
% QM = load('QDSBN1_18152');    % 18152 Note8 B-1
% QM = load('QDSBN2_18152');    % 18152 Note8 B-2
% QM = load('QDSBN3_18152');    % 18152 Note8 B-3
% QM = load('QDSBN4_18152');    % 18152 Note8 B-4
% QM = load('QDSBN_18156');    % 18156 Note8
% QM = load('QDSBN1_18156');    % 18156 Note8 B-1
% QM = load('QDSBN2_18156');    % 18156 Note8 B-2
% QM = load('QDSBN3_18156');    % 18156 Note8 B-3
% QM = load('QDSBN4_18156');    % 18156 Note8 B-4
% QM = load('QDSBS_18156');    % 18156 S8
% QM = load('QDSBS1_18156');    % 18156 S8 B-1
% QM = load('QDSBS2_18156');    % 18156 S8 B-2
% QM = load('QDSBS3_18156');    % 18156 S8 B-3
% QM = load('QDSBS4_18156');    % 18156 S8 B-4
% QM = load('QIHUN_18157');    % 18157 Note8 Inha 30분
% QM = load('QIHUN_18157_whole');    % 18157 Note8 Inha whole
% QM = load('QDSBN1_18157');    % 18157 N8 B-1
% QM = load('QDSBN2_18157');    % 18157 N8 B-2
% QM = load('QDSBS1_18157');    % 18157 S8 B-1
% QM = load('QDSBS2_18157');    % 18157 S8 B-2
% QM = load('QDSBS_18163_geo_3');    % 18163 S8 Geo
% QM = load('QDSBN2_18163');    % 18163 N8 gnss logger
% QM = load('QDSBS9_18164');    % 18164 S9 gnss logger
% QM = load('QDSBS_18164_1');    % 18164 S9 gnss logger(교수님) 1
% QM = load('QDSBS_18164_2');    % 18164 S9 gnss logger(교수님) 2
% QM = load('QDBS9_18165_b');    % 18165 S9 gnss logger b
% QM = load('QDBS9_18165_a');    % 18165 S9 gnss logger b
load('18169_rt_data.mat');
% QM = QDAN8_18169;    % 18165 JAVAD-A
QM = QDBUB_18169;    % 18165 S9(GEO++) L1 : cycle
QM(:,1) = round(QM(:,1));
wholeTTs = unique(QM(:,1));
%% 존재 위성군 판단
Constell = [length(find(QM(:,3) < 200)), length(find(QM(:,3) > 200 & QM(:,3) < 300)),...
    length(find(QM(:,3) > 300 & QM(:,3) < 400)), length(find(QM(:,3) > 400 & QM(:,3) < 500)),...
    length(find(QM(:,3) > 500 & QM(:,3) < 600)), length(find(QM(:,3) > 700))];
ExistCons = find(Constell(:) ~= 0);

%% 시계열 전체 수집 위성수
for i = 1:length(wholeTTs)
    gs = wholeTTs(i);
    numsats = QM(find(QM(:,1) == gs), 2);
    numsats = length(unique(numsats));
    whole(i,:) = [gs, numsats];
end
tHour = mod(whole(:,1), 86400); tHour = tHour/3600;
if find(tHour(:) == 0) > 1
    tHour(find(tHour(:) == 0):end) = tHour(find(tHour(:) == 0):end) + 24;
end
figure(1)
stairs(tHour, whole(:,2) ,'bo-')
hold on;


for C=1:length(ExistCons)
    Cons = ExistCons(C);
    Type = Cons * 100;
    QM_C1 = QM(QM(:,3) == Type + 20, :);        % Base C1 데이터 추출
    QM_L1 = QM(QM(:,3) == Type + 11, :);        % Base L1 데이터 추출
    prnlist = unique(QM_C1(:,2));
    FinalTTs = unique(QM(:,1));
    tHour = mod(FinalTTs, 86400); tHour = tHour/3600;
    if find(tHour(:) == 0) > 1
        tHour(find(tHour(:) == 0):end) = tHour(find(tHour(:) == 0):end) + 24;
    end
    
    Culoum = zeros(length(prnlist),200);
    Culoum_gs = zeros(length(prnlist),200);
    Culoum_diff = zeros(length(prnlist),200);
    %% 위성군/위성별 L1 시계열 plot
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
    %% 위성군 별 위성수 시계열 plot
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
    
    figure(1)
    stairs(numprn(:,1), numprn(:,2),'o-')
    grid on; hold on;
    xlim([min(numprn(:,1)), max(numprn(:,1))]);
    xlabel('Hour')
    figure(Type)
    stairs(numprn(:,1), numprn(:,2),'bo-')
    grid on;
    xlim([min(numprn(:,1)), max(numprn(:,1))]);
    xlabel('Hour')
    ylim([1, 12])
    mean(numprn(:,2))
    numprn = [];
end

figure(1)
legend('Total','GPS','BDS','GLONASS','QZSS')


% % figure(199)
% % stairs(states(:,1), states(:,2),'bo-')
% % grid on;
% % xlim([min(states(:,1)), max(states(:,1))]);
% % xlabel('Hour')
% % ylim([0, 5])