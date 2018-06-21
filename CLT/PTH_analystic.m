clear all; close all

%% load
load('test1_L.mat');
test1 = estm(find(estm(:,1) > 105000 & estm(:,1) < 105343),:);
%% 이상현상 발생 QMfile Load
FileQM='test1_L';
[arrQM_1_L, FinalPRNs_1_L, FinalTTs_1_L] = ReadQM(FileQM);
% arrQM_1_L = arrQM_1_L(find(arrQM_1_L(:,1) > 1000),:);
arrQM_1_L = arrQM_1_L(find(arrQM_1_L(:,1) > 105000 & arrQM_1_L(:,1) < 105343),:);
FileQM='test1_R';
[arrQM_1_R, FinalPRNs_1_R, FinalTTs_1_R] = ReadQM(FileQM);
% arrQM_1_R = arrQM_1_R(find(arrQM_1_R(:,1) > 1000),:);
arrQM_1_R = arrQM_1_R(find(arrQM_1_R(:,1) > 105000 & arrQM_1_R(:,1) < 105343),:);
% FileQM='test2_1_L';
% [arrQM_2_1_L, FinalPRNs_2_1_L, FinalTTs_2_1_L] = ReadQM(FileQM);
% FileQM='test2_2_R';
% [arrQM_2_2_R, FinalPRNs_2_2_R, FinalTTs_2_2_R] = ReadQM(FileQM);
% FileQM='test4_1_L';
% [arrQM_4_1_L, FinalPRNs_4_1_L, FinalTTs_4_1_L] = ReadQM(FileQM);
% 
FinalPRNs_1_L_GPS = unique(arrQM_1_L(find(arrQM_1_L(:,3) ==120),2));
FinalPRNs_1_L_BDS = unique(arrQM_1_L(find(arrQM_1_L(:,3) ==220),2));
% 
tHour = mod(arrQM_1_L(:,1),86400)/3600;
arrQM_1_L(:,1) = tHour;
tHour = mod(arrQM_1_R(:,1),86400)/3600;
arrQM_1_R(:,1) = tHour;
for i = 1:length(FinalPRNs_1_L_GPS)
    prn_gps = FinalPRNs_1_L_GPS(i);
    prn_gps_pr = arrQM_1_L(find(arrQM_1_L(:,2) == prn_gps & arrQM_1_L(:,3) == 120),:);
    prn_gps_ = prn_gps_pr(2:end,:);
    Prn_gps = prn_gps_pr(1:end-1,:) - prn_gps_ ;
    prn_gps_snr = arrQM_1_L(find(arrQM_1_L(:,2) == prn_gps & arrQM_1_L(:,3) == 141),:);
    prn_gps_snr_ = prn_gps_snr(2:end,:);
    Prn_gps_snr = prn_gps_snr(1:end-1,:) - prn_gps_snr_ ;
%     prn_gps_pr_R = arrQM_1_R(find(arrQM_1_R(:,2) == prn_gps & arrQM_1_R(:,3) == 120),:);
%     prn_gps_snr_R = arrQM_1_R(find(arrQM_1_R(:,2) == prn_gps & arrQM_1_R(:,3) == 141),:);
%     figure(prn_gps)
figure(prn_gps)
suptitle(['prn = ', num2str(prn_gps)])
    
    subplot(3,1,1)
    title('X(ecef)')
    plot(mod(test1(:,1),86400)/3600,test1(:,2),'b-')
%     plot(prn_gps_snr(:,1), prn_gps_snr(:,4),'b:')
    hold on; grid on;
    xlim([min(arrQM_1_L(:,1)), max(arrQM_1_L(:,1))])

    subplot(3,1,2)
    title('\Delta Pseudo-range')
    plot(prn_gps_pr(:,1), prn_gps_pr(:,4),'b.')
%     plot(prn_gps_(:,1)-0.5, Prn_gps(:,4),'b-')
    hold on; grid on;
%     plot(prn_gps_pr_R(:,1), prn_gps_pr_R(:,4),'ro')
%     xlim([min(tHour), max(tHour)])
    xlim([min(arrQM_1_L(:,1)), max(arrQM_1_L(:,1))])
    
    subplot(3,1,3)
    title('\Delta SNR')
    plot(prn_gps_snr(:,1), prn_gps_snr(:,4),'b-')
%     plot(prn_gps_snr_(:,1), Prn_gps_snr(:,4),'b-')
    hold on; grid on;
    xlim([min(arrQM_1_L(:,1)), max(arrQM_1_L(:,1))])
    
%     subplot(4,1,4)
%     title('Number of Sats')
%     plot(test1(:,1),test1(:,7),'b-')
%     hold on; grid on;
%     ylim([4 10])
%     xlim([min(arrQM_1_L(:,1)), max(arrQM_1_L(:,1))])
end
% 
% % for i = 1:length(FinalPRNs_1_L_BDS)
% %     prn_bds = FinalPRNs_1_L_BDS(i);
% %     prn_bds_pr = arrQM_1_L(find(arrQM_1_L(:,2) == prn_bds & arrQM_1_L(:,3) == 220),:);
% %     prn_bds_snr = arrQM_1_L(find(arrQM_1_L(:,2) == prn_bds & arrQM_1_L(:,3) == 241),:);
% %     figure(prn_bds)
% %     subplot(2,1,1)
% %     plot(prn_bds_pr(:,1), prn_bds_pr(:,4),'bo')
% %     xlim([min(tHour), max(tHour)])
% % %     xlim([min(arrQM_1_L(:,1)), max(arrQM_1_L(:,1))])
% %     subplot(2,1,2)
% %     plot(prn_bds_snr(:,1), prn_bds_snr(:,4),'b-')
% % %     xlim([min(arrQM_1_L(:,1)), max(arrQM_1_L(:,1))])
% %     xlim([min(tHour), max(tHour)])
% % end