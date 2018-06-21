%
%   coded by Joonseong Gim, Jan 13, 2016
%

clear all; close all;


%% 불변 변수 설정: 빛의 속도, 관측치
CCC = 299792458.;   % CCC = Speed of Light [m/s]
ObsType = 120;      % 사용할 관측치 설정 - 120: C/A = C1
ObsType2 = 141;      % 사용할 관측치 설정 - 141: snr = S1
TruePos = [-3041235.578 4053941.677 3859881.013];   % : JPspace A point
% TruePos = [-3041241.741 4053944.143 3859873.640];   % : JPspace B point
%% 임계고도각 설정
eleCut = 15;

%% rinex file 입력
% obsfile = 'DDBs115_37.16o';
% navfile = 'DDBs115_37.16n';
% obsfile = '20160211A_2.obs';
% navfile = '20160211A_2.nav';
% obsfile = 'r1.16o';
% navfile = 'r1.16n';
% obsfile = 'l1.obs';
% navfile = 'l1.nav';
obsfile = 'SDT1_160824.obs';
navfile = 'brdc2370.16n';
% obsfile = 'DAUU049r.16.obs';
% navfile = 'brdc0490.16n';
%% Observation Rinex 파일로 부터 QMfile 생성
WriteObs(obsfile)

%% 만들어진 QMfile을 obsevation Rinex 파일 관측 날짜, 시간으로 변경
% rename = renameQMfile(obsfile);
% rename = 'QMDBUU049i_16';
rename = 'QSDT1_16_ob082';
[YY, DOY] = obs2YYDOY(obsfile);
[gw, GD] = ydoy2gwgd(YY, DOY); %: GPS WEEK 결정
%% QM 파일 읽어들여서 행렬로 저장하고, 사용할 관측치 추출
[arrQM, FinalPRNs, FinalTTs] = ReadQM(rename);
QM = SelectQM(arrQM, ObsType);
QM2 = SelectQM(arrQM, ObsType2);

%% 항법메시지를 읽어들여서 행렬로 저장하고, Klobuchar 모델 추출
eph = ReadEPH(navfile);
[al, be] = GetALBE(navfile);
%% 라이넥스 파일에서 대략적인 관측소 좌표를 뽑아내고 위경도로 변환
% AppPos = [-3041241.741 4053944.143 3859873.640];   % : JPspace B point
AppPos = GetAppPos(obsfile);
if AppPos(1) == 0
    AppPos = TruePos;
end

gd = xyz2gd(AppPos); AppLat = gd(1); AppLon = gd(2);

%% 추정에 필요한 초기치 설정
MaxIter = 4;
EpsStop = 1e-5;
ctr = 1;
x = [AppPos ctr]; x = x';

%% 추정과정 시작
NoEpochs = length(FinalTTs);
nEst = 0;


for j = 1:NoEpochs
    indexQM = find(QM(:,1) == FinalTTs(j));
    QM_1 = QM(indexQM,:);
    QM_2 = QM2(indexQM,:);
    existprn = intersect(unique(eph(:,18)), QM_1(:,2));
    arrSV = zeros(length(existprn),1);
    for k = 1:length(existprn)
        arrSV(k) = find(QM_1(:,2) == existprn(k,1));
    end
    QM_1 = QM_1(sort(arrSV),:);
    QM_2 = QM_2(sort(arrSV),:);
    if length(existprn) < 5
        nEst = nEst + 1;
        estm(nEst,:) = [0,0,0,0,0];
    else
        for Iter = 1:MaxIter
            HTH = zeros(4,4);
            HTy = zeros(4,1);
            
            NoSats = length(QM_1);
            gs = QM_1(1,1);
            
            vec_site = x(1:3)';
            visiSat(j,1) = gs; visiSat(j,2) = NoSats;
            ZHD = TropGPTh(vec_site, gw, gs);                 %: TROP: GPT
            
            for i = 1:NoSats
                prn = QM_1(i,2);
                obs = QM_1(i,4);
                S1 = QM_2(i,4);
                icol = PickEPH(eph, prn, gs);
                toe = eph(icol, 8); a = eph(icol, 19); b = eph(icol, 20); c = eph(icol, 21); Tgd = eph(icol, 23);
                %----- 신호전달시간 계산
                STT = GetSTTbrdc(gs, prn, eph, vec_site);
                tc = gs - STT;
                %----- 위성궤도 계산
                vec_sat = GetSatPosNC(eph, icol, tc);
                vec_sat = RotSatPos(vec_sat, STT);                      %: 지구자전 고려
                %----- 최종 RHO 벡터 계산
                vec_rho = vec_sat - vec_site;
                rho = norm(vec_rho);
                [az,el] = xyz2azel(vec_rho, AppLat, AppLon);
                
                if el >= eleCut %15
                    %                     W = MakeW_elsnr(el,S1);
                    %                     W = MakeW_elpr(el);
                    W = 1;
                    dRel = GetRelBRDC(eph, icol, tc);
                    dtSat = a + b*(tc - toe) + c*(tc - toe)^2 - Tgd + dRel;
                    dIono = ionoKlob(al, be, gs, az, el, vec_site);
                    %                 dTrop_H = Hopfield(el, 11, vec_site, 9999);                   % Hopfield Model
                    %                 dTrop_S = Saastamoinen(el, 11, vec_site, 9999);               % Saastamoinen Model
                    dTrop_G = ZHD2SHD(gw, gs, vec_site, el, ZHD);                   % GPT model
                    %                 com = rho + x(4) - CCC * dtSat + dIono + dTrop_H;             % Hopfield Model
                    %                 com = rho + x(4) - CCC * dtSat + dIono + dTrop_S;             % Saastamoinen Model
                    com = rho + x(4) - CCC * dtSat + dIono + dTrop_G;             % GPT Model
                    %                     com = rho + x(4) - CCC * dtSat;
                    y = obs - com;
                    H = [ -vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 1];
                    HTH = HTH + H'*W*H;
                    HTy = HTy + H'*W*y;
                    
                end
            end
            xhat = inv(HTH) * HTy;
            x = x + xhat;
            
            if norm(xhat) < EpsStop;
                nEst = nEst + 1;
                estm(nEst,1) = gs;
                estm(nEst,2:5) = x(1:4);
                fprintf('gs: %6.0f     %2.5f \n',gs,norm(TruePos'-x(1:3)));
                break;
            end
        end
    end
    user_gd(j,:) = xyz2gd(estm(j,2:4)); % user의 xyz값 gd 로 변환
    %         AppLat = user_gd(j,1); AppLon = user_gd(j,2);
    %         user_xyz(j,:) = estm(j,2:4);        % user xyz값 행렬로 변환
    %
end
estm = estm(find(estm(:,1) > 0),:);
for es = 1:length(estm(:,1))
    user_gd(es,:) = xyz2gd(estm(es,2:4)); % user의 xyz값 gd 로 변환
    AppLat = user_gd(es,1); AppLon = user_gd(es,2);
    user_xyz(es,:) = estm(es,2:4);        % user xyz값 행렬로 변환
end
%% user_position mean값 기준 topology 계산 과정
% [user_mean] = PlotMeanTopo(user_xyz);
% user_mean = [mean(user_xyz(:,1)) mean(user_xyz(:,2)) mean(user_xyz(:,3))];
% user_mean = [-3055625.93798130,4034964.65967131,3868139.09930058];
% user_mean = [-3055530.33388512,4035091.49810586,3868085.81846658];
% user_mean = TruePos;
% [dXYZ, dNEV] = PosTErrorsJOON(estm(:, 1), user_mean, estm(:, 2:5));

% [dXYZ, dNEV] = PosTErrorsJOON(estm(:, 1), TruePos, estm(:, 2:5));
% [dXYZ, dNEV] = PosTErrors(estm(:, 1), TruePos, estm(:, 2:5));

%% QM type 별 plot
% % PPPlotQM(rename,141)


%% 구글 plot
% figure(201)
% axis([127.03 127.041 37.528 37.548]);
% plot_google_map;
% axis equal
% axis([127.03 127.041 37.528 37.548]);
%
% for i = 1:length(user_gd)
%     GD_lon(i,1) = user_gd(i,2);
%     GD_la(i,1) = user_gd(i,1);
%     figure(201)
%     grid on
%     hold on
%     plot(GD_lon, GD_la,'ro','markeredgecolor','y','markerfacecolor','r','markersize',3)
%     drawnow
% end

VRS_text = 'SDT1_VRS_16237.txt';
[target,UBLOX] = gapconv(VRS_text, 0.43, 16, 237);

%% PostErrors
[dXYZ, dNEV] = PosTErrors5(estm, UBLOX, visiSat);
% 
%% 위성별 SNR
% PPPlotQM(QM2,141)

% figure(201)
% grid on
% hold on
% % axis([(min(user_gd(:,2))-0.0001) (max(user_gd(:,2))+0.0001) (min(user_gd(:,1))-0.0001) (max(user_gd(:,1))+0.0001)]);
% plot(user_gd(:,2), user_gd(:,1),'ro','markeredgecolor','y','markerfacecolor','r','markersize',3)
% plot_google_map;