clear all
close all

% States = load('QDSBS_18135_state');      % 18135 S8
% States = load('QDSBS_18137_state');      % 18137 S8
% States = load('QDSBS_18145_state');      % 18145 S8
% States = load('QDSBS_18147_state');      % 18147 S8
% States = load('QDSBS1_18148_state');      % 18148 S8 1
% States = load('QDSBS2_18148_state');      % 18148 S8 2
% States = load('QDSBN2_18148_state');      % 18148 N9 2
% States = load('QIHU4_18149_state');    % 18149 INHA Note8
% States = load('QDSBNU_18149_state');    % 18149 NEXUS9 up
% States = load('QDSBND_18149_state');    % 18149 NEXUS9 down
States = load('QDSBS9_18164_state');      % 18148 S8 2
prnlist = unique(States(:,2));
FinalTTs = unique(States(:,1));
% FinalTTs = round(unique(States(:,1)));
% States(:,1) = round(States(:,1));
tHour = mod(FinalTTs, 86400); tHour = tHour/3600;
if find(tHour(:) == 0) > 1
    tHour(find(tHour(:) == 0):end) = tHour(find(tHour(:) == 0):end) + 24;
end

for i=1:length(prnlist)
prn = prnlist(i);
states = States(find(States(:,2) == prn),:);
prn_tHour = mod(states(:,1), 86400); prn_tHour = prn_tHour/3600;
    if find(prn_tHour(:) == 0) > 1
        prn_tHour(find(prn_tHour(:) == 0):end) = prn_tHour(find(prn_tHour(:) == 0):end) + 24;
    end
figure(299)
plot(prn_tHour, states(:,3))
hold on; grid on;
xlim([min(tHour), max(tHour)])
ylim([0,5])
figure(prn)
plot(prn_tHour, states(:,3))
hold on; grid on;
xlim([min(tHour), max(tHour)])
ylim([0,5])
end

