clear all
close all

Dynamic = 0;
SYS = 4;
%% VRS load
% Base_vrs = load('DD1Bs_adm.txt');
% Rover_vrs = load('DD1Rv_adm.txt');
% Base_vrs = load('DD2Bs_adm.txt');
% Rover_vrs = load('DD2Rv_adm.txt');
% Base_vrs = load('PTCO4_joon_170216_adm.txt');
% Rover_vrs = load('PTCO4_hyunu_170216_adm.txt');
%% 좌표 참값
Truedis = 9.9209;
% Truedis = 1.1;
% TruePos1 = [-3041235.57800000,4053941.67700000,3859881.01300000];        % 대성 A
% TruePos2 = [-3041241.74100000,4053944.14300000,3859873.64000000];        % 대성 B
TruePos = [-3041235.57800000,4053941.67700000,3859881.01300000];        % 대성 A
% TruePos = [-3041241.74100000,4053944.14300000,3859873.64000000];        % 대성 B
% TruePos = [-3041092.515 4054082.611 3859643.389];        % F point
% TruePos = [-3027386.463213997 4071581.638074351 3852036.292033684]; % 현우집앞
% TruePos = [-3108706.97103816,4078522.84147406,3779757.23816543];    % 오송 1
% TruePos = [-3108697.15866998,4078501.37490046,3779789.12574991];    % 오송 2
% TruePos = [-3051683.19233025,4044464.46325792,3861660.30965950];    % 강남 2
% TruePos = [-3080121.50057624,4058693.37773498,3824124.32021532];    % 안성W 1-1
% obsfile = '170509_r2.obs';
% TruePos = [-3053365.29481677,4039290.16703344,3865445.80715444];      % 집앞
% TruePos = App_pos(obsfile);
%% eph load
% load('eph170125.mat'); DOY = 025; YY  = 17;        % 대성
% load('eph170216.mat'); DOY = 047; YY  = 17;        % 신항대로
% load('eph170322.mat'); DOY = 081; YY  = 17;        % 오송
% load('eph170323.mat'); DOY = 082; YY  = 17;        % 강남
% load('eph170508.mat'); DOY = 128; YY  = 17;        % 강남2
% load('eph170802.mat'); DOY = 214; YY  = 17;        % 안성W
load('eph180228.mat'); DOY = 059; YY  = 18;        % 2018-02-28 teheran
% load('eph180305.mat'); DOY = 064; YY  = 18;        % 2018-03-05 F point
load('eph180618.mat');
navfile = strcat('brdm', num2str(DOY,'%03d'), '0.', num2str(YY,'%02d'), 'p');   %: Navigation RINEX file
gps_nav = strcat(navfile(1:3),'c',navfile(5:8),'.',navfile(10:11),'n');
[al, be] = GetALBE(gps_nav);

[gw, gd] = ydoy2gwgd(YY, DOY); %: GPS WEEK 결정

%% QMfile load
% QMfileBs = 'QM170125_A';        % 대성 A
% QMfileRv = 'QM170125_B';        % 대성 B
% QMfileBs = 'ublox_joon';
% QMfileRv = 'ublox_hyunu';
% QMfileBs = 'QM170322_Bs_1';        % 오송 bs 1
% QMfileRv = 'QM170322_Rv_1';        % 오송 rv 1
% QMfileBs = 'QM170322_Bs_2';        % 오송 bs 2
% QMfileRv = 'QM170322_Rv_2';        % 오송 rv 2
% QMfileBs = 'QM170323_Bs';        % 강남 bs 
% QMfileRv = 'QM170323_Rv';        % 강남 rv 
% QMfileBs = 'QM170509_Bs3';        % 강남 bs 
% QMfileRv = 'QM170509_Rv3';        % 강남 rv 
% QMfileBs = 'ANSG_17214_u1_1'; QMfileRv = 'ANSG_17214_u1_2'; 
% QMfileBs = 'QTHU1_18059';        % 2018-02-28 teheran u-blox1
QMfileRv = 'QTHU2_18059';        % 2018-02-28 teheran u-blox2
QMfileRv = 'QDBUB_18169';        % 2018-02-28 teheran u-blox2
% QMfileBs = 'QTHJ1_18059';        % 2018-02-28 teheran javad1
% QMfileRv = 'QTHJ2_18059';        % 2018-02-28 teheran javad2
% QMfileBs = 'QTHS1_18059';        % 2018-02-28 teheran septentrio1
% QMfileRv = 'QTHS2_18059';        % 2018-02-28 teheran septentrio2
% QMfileBs = 'QS8U1_18064_';        % 2018-03-05 F point S8
% 
% Base_PP = PP_g(QMfileBs,eph,TruePos,DOY,YY);              % with Correction
estm_PP = PP_g(QMfileRv,eph,TruePos,DOY,YY);              % with Correction
% Base_PP= PP_g_woc(QMfileBs,eph,TruePos,DOY,YY);              % with Correction
% estm_PP= PP_g_woc(QMfileRv,eph,TruePos,DOY,YY);              % with Correction
% Base_PP = PP_c(QMfileBs,eph,TruePos,DOY,YY);              % with Correction
% estm_PP = PP_c(QMfileRv,eph,TruePos,DOY,YY);              % with Correction
% Base_PP= PP_c_woc(QMfileBs,eph,TruePos,DOY,YY);              % with Correction
% estm_PP= PP_c_woc(QMfileRv,eph,TruePos,DOY,YY);              % with Correction
% Base_PP= PP_r(QMfileBs,eph,TruePos,DOY,YY);              % with Correction
% estm_PP= PP_r(QMfileRv,eph,TruePos,DOY,YY);              % with Correction
% Base_PP= PP_r_woc(QMfileBs,eph,TruePos,DOY,YY);              % with Correction
% estm_PP= PP_r_woc(QMfileRv,eph,TruePos,DOY,YY);              % with Correction
% Base_PP= PP_gc_old(QMfileBs,eph,TruePos,DOY,YY);              % with Correction
% estm_PP= PP_gc_old(QMfileRv,eph,TruePos,DOY,YY);              % with Correction
% Base_PP= PP_gc_dop(QMfileBs,eph,TruePos,DOY,YY);              % with Correction
% estm_PP= PP_gc_dop(QMfileRv,eph,TruePos,DOY,YY);              % with Correction
% Base_PP= PP_gc_kf2(QMfileBs,eph,TruePos1,DOY,YY);              % with Correction
% estm_PP= PP_gc_kf2(QMfileRv,eph,TruePos2,DOY,YY);              % with Correction
% Base_PP= PP_gc_kf2(QMfileBs,eph,TruePos,DOY,YY);              % with Correction
% estm_PP= PP_gc_kf2(QMfileRv,eph,TruePos,DOY,YY);              % with Correction
% Base_PP= PP_gc_woc(QMfileBs,eph,TruePos,DOY,YY);              % with Correction
% estm_PP= PP_gc_woc(QMfileRv,eph,TruePos,DOY,YY);              % with Correction
% Base_PP= PP_gr(QMfileBs,eph,TruePos,DOY,YY);              % with Correction
% estm_PP= PP_gr(QMfileRv,eph,TruePos,DOY,YY);              % with Correction
% Base_PP= PP_gr_woc(QMfileBs,eph,TruePos,DOY,YY);              % with Correction
% estm_PP= PP_gr_woc(QMfileRv,eph,TruePos,DOY,YY);              % with Correction
% Base_PP= PP_rc(QMfileBs,eph,TruePos,DOY,YY);              % with Correction
% estm_PP= PP_rc(QMfileRv,eph,TruePos,DOY,YY);              % with Correction
% Base_PP= PP_rc_woc(QMfileBs,eph,TruePos,DOY,YY);              % with Correction
% estm_PP= PP_rc_woc(QMfileRv,eph,TruePos,DOY,YY);              % with Correction
% Base_PP= PP_grc(QMfileBs,eph,TruePos,DOY,YY);              % with Correction
% estm_PP= PP_grc(QMfileRv,eph,TruePos,DOY,YY);              % with Correction
% Base_PP= PP_grc_woc(QMfileBs,eph,TruePos,DOY,YY);              % with Correction
% estm_PP= PP_grc_woc(QMfileRv,eph,TruePos,DOY,YY);              % with Correction

% estm = estm(650:end, :);
% Base = Base(650:end, :);

% if DOY == 47 & YY == 17
%     estm_PP = estm_PP(find(estm_PP(:,1) > 387842),:);
%     Base_PP = Base_PP(find(Base_PP(:,1) > 387842),:);
% end
% %% dynamic에 따라 그래프 그리는 방법 선택
% if Dynamic == 0
%     [DDdXYZ, DDdNEV, DDdis, DDrms, DDstd] = PostErrorsDDstatic(estm_PP, Base_PP, Truedis,8, 12, SYS);    % 임의의 장소에서 이동 측위시
% elseif Dynamic == 1
%     [DDdXYZ, DDdXYZ_vrs, DDdNEV, DDdNEV_vrs, DDdis, DDrms, result] = PostErrorsDDkine(estm_PP, Base_PP, Base_vrs, Rover_vrs, -2, 2, SYS);    % 임의의 장소에서 이동 측위시
% end

%% 대성 AB일때 그래프 선택
% if DOY == 25 & YY == 17
%     [dXYZ, dNEV, result] = PostErrorskine(estm_PP(800:end,:), [-3041241.74100000,4053944.14300000,3859873.64000000], SYS); 
% %     [dXYZ, dNEV, result] = PostErrorskine(estm_PP(800:4400,:), [-3041241.74100000,4053944.14300000,3859873.64000000], SYS); 
% %     [dXYZ, dNEV, result] = PostErrorskine(estm_PP(800:2600,:), [-3041241.74100000,4053944.14300000,3859873.64000000], SYS); 
% %     [dXYZ, dNEV, result] = PostErrorskine(estm_PP(2601:4400,:), [-3041241.74100000,4053944.14300000,3859873.64000000], SYS); 
% %     [dXYZ, dNEV, result] = PostErrorskine(estm_PP(2401:3600,:), [-3041241.74100000,4053944.14300000,3859873.64000000], SYS); 
%     
%     [dXYZ, dNEV, result] = PostErrorskine(Base_PP(800:end,:), [-3041235.57800000,4053941.67700000,3859881.01300000], SYS); 
% %     [dXYZ, dNEV, result] = PostErrorskine(Base_PP(800:4400,:), [-3041235.57800000,4053941.67700000,3859881.01300000], SYS); 
% %     [dXYZ, dNEV, result] = PostErrorskine(Base_PP(800:2600,:), [-3041235.57800000,4053941.67700000,3859881.01300000], SYS); 
% %     [dXYZ, dNEV, result] = PostErrorskine(Base_PP(2601:4400,:), [-3041235.57800000,4053941.67700000,3859881.01300000], SYS); 
% %     [dXYZ, dNEV, result] = PostErrorskine(Base_PP(2401:3600,:), [-3041235.57800000,4053941.67700000,3859881.01300000], SYS);
%     
% elseif DOY == 82
% else
%     [dXYZ, dNEV, result] = PostErrorskine(estm_PP, Rover_vrs, SYS);
% end

%% 구글 plot
% Base_PP(:,1) = round(Base_PP(:,1));
% estm_PP(:,1) = round(estm_PP(:,1));
% FinalTTs = intersect(Base_PP(:,1), estm_PP(:,1));
% for i = 1:length(FinalTTs)
%     gs = FinalTTs(i); 
%     Base_gd(i,1:4) = [gs xyz2gd(Base_PP(find(Base_PP(:,1) == gs),2:4))];
%     estm_gd(i,1:4) = [gs xyz2gd(estm_PP(find(estm_PP(:,1) == gs),2:4))];
% end

% % base plot
% FinalTTs = Base_PP(:,1);
% for i = 1:length(FinalTTs)
%     gs = FinalTTs(i); 
%     Base_gd(i,1:4) = [gs xyz2gd(Base_PP(find(Base_PP(:,1) == gs),2:4))];
% %     estm_gd(i,1:4) = [gs xyz2gd(estm_PP(find(estm_PP(:,1) == gs),2:4))];
% end
% figure(10)
% hold on; grid on;
% plot(Base_gd(:,3), Base_gd(:,2), 'b.-')
% % xlim([127.00365234782228, 127.072695849665])
% % ylim([37.479551317952414,37.524398954977315])
% plot_google_map

% rover plot
FinalTTs = estm_PP(:,1);
for i = 1:length(FinalTTs)
    gs = FinalTTs(i); 
%     Base_gd(i,1:4) = [gs xyz2gd(Base_PP(find(Base_PP(:,1) == gs),2:4))];
    estm_gd(i,1:4) = [gs xyz2gd(estm_PP(find(estm_PP(:,1) == gs),2:4))];
end
figure(20)
hold on; grid on;
plot(estm_gd(:,3), estm_gd(:,2), 'r.')
plot_google_map

[dXYZ, dNEV] = PosTErrorsJOON(Base_PP(:,1), TruePos, Base_PP(:,2:4))
% figure()
% hold on; grid on;
% plot(Base_gd(:,3), Base_gd(:,2), 'bo')
% plot(estm_gd(:,3), estm_gd(:,2), 'ro')
% plot_google_map
% mfilename