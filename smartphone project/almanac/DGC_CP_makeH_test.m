% %
% %   modified by Joonseong Gim, aug 7, 2016
% %
% 
% clear all; close all;
% 
% %% QMfile 호출
% FileQM = 'QM170125_A';
% FileQM = 'QM170811_B';
% 
% %% Almanac file & almanac matrix
% load('almanac_17255.mat');
% alm = [alm_gps;alm_bds];
% 
% %% eph load
% % load('eph170125.mat'); DOY = 025; YY  = 17;        % 대성
% load('eph170811.mat'); DOY = 223; YY  = 17;        % 대성
% 
% %% load a PRC file
% PRCfile = 'PPS1_170125.t41';
% PRCfile = 'PPS1_170811.t41';
% RTCM = load('PPS1_170125.t41');
% %% 불변 변수 설정: 빛의 속도, 관측치
% CCC = 299792458.;   % CCC = Speed of Light [m/s]
% %% 임계고도각 설정
% eleCut = 15;
% %% 윤초 핸들링
% LeapSecBDS = 14;
% 
% %% 사용할 관측치 선정
% g_ObsType = 120; % gps C1
% g_ObsType_snr = 141;
% c_ObsType = 220; % bds C1
% c_ObsType_snr = 241;
% 
% %% 기준좌표 설정
% % TruePos = [-3041235.578 4053941.677 3859881.013];   % : JPspace A point
% TruePos = [-3041241.741 4053944.143 3859873.640];   % : JPspace B point
% 
% %% 만들어진 QMfile을 obsevation Rinex 파일 관측 날짜, 시간으로 변경
% % rename = renameQMfile(obsfile);
% [gw, GD] = ydoy2gwgd(YY, DOY); %: GPS WEEK 결정
% 
% %% QM 파일 읽어들여서 행렬로 저장하고, 사용할 관측치 추출
% [arrQM, FinalPRNs, FinalTTs] = ReadQM(FileQM);
% arrQM = arrQM(find(arrQM(:,3) < 300),:);
% arrQM(:,1) = round(arrQM(:,1));
% 
% QM = SelectQM_gc(arrQM, g_ObsType, c_ObsType);
% QM_snr = SelectQM_gc(arrQM, g_ObsType_snr, c_ObsType_snr);
% FinalTTs = unique(QM(:,1));
% 
% [GPSPRC, BDSPRC, GLOPRC, PRC_sorted] = PRCsort(PRCfile, QM);
% 
% %% 항법메시지를 읽어들여서 행렬로 저장하고, Klobuchar 모델 추출
% gps_nav = strcat('brdc', num2str(DOY,'%03d'), '0.', num2str(YY,'%02d'), 'n');   %: Navigation RINEX file
% fid = fopen(gps_nav,'r');
% if fid == -1
%     al = zeros(4,1); be = zeros(4,1);
% else
%     [al, be] = GetALBE(gps_nav);
% end
% 
% AppPos = TruePos;
% 
% gd = xyz2gd(AppPos); AppLat = gd(1); AppLon = gd(2);
% 
% %% 추정에 필요한 초기치 설정
% Maxiter = 10;
% EpsStop = 1e-5;
% ctr = 1; ctr2 = 1;
% x = [AppPos ctr ctr2]; x = x';
% x_cp = [AppPos ctr ctr2]; x_cp = x_cp';
% 
% %% 추정과정 시작
% NoEpochs = length(FinalTTs);
% EstPos = zeros(NoEpochs,5);
% nEst = 0;
% j=1;
% estm = zeros(NoEpochs,6);
% 
% cnt = 1;

clear all
close all


% load('DGNSS_CP_alm_test170125.mat');
load('DGNSS_CP_alm_test170811.mat');
x_cp_alm = x;
x_dgc = x;
tic
load('eph1701252.mat'); 
% for j = 1:NoEpochs
for j = 1:600
    

    gs = FinalTTs(j);
    %     QMe_ = qmHandle(QMd.pickQM(gs,':',':'));
    %     QMe = dopplerSM(QMe_.getQM(),5);
    indexQM = find(QM(:,1) == gs);
    QM1e = QM(indexQM,:);
    QM1e_snr = QM_snr(indexQM,:);
    existprn = intersect(QM1e(:,2), unique(eph(find(eph(:,22) == 0),18)));
    NoSats = length(existprn);
    
    for iter = 1:Maxiter
        
        if NoSats <= 6
            break
        end
        vec_site = x(1:3)';
        AppLLH = xyz2gd(vec_site); AppLat = AppLLH(1); AppLon = AppLLH(2);
        vec_site_dgc = x_dgc(1:3)';
        vec_site_cp = x_cp(1:3)';
        
        ZHD = TropGPTh(TruePos, gw, gs); %: TROP: GPT
        NoUsed_gps = 0;
        NoUsed_bds = 0;
        NoGPSsUsed = 0;
        NoBDSsUsed = 0;
        NoGLOsUsed = 0;
        cnt2=1;
        
        for i = 1:NoSats
            prn = QM1e(find(QM1e(:,2) == existprn(i)),2);
            obs = QM1e(find(QM1e(:,2) == existprn(i)),4);
            snr = QM1e_snr(find(QM1e_snr(:,2) == existprn(i)),4);
            icol = PickEPH(eph, prn, gs);
            
            toe = eph(icol, 8); a = eph(icol, 19); b = eph(icol, 20); c = eph(icol, 21); Tgd = eph(icol, 23);
            %             STT=0;
            STT = GetSTTbrdm(gs, eph, icol, x(1:3)); % 신호전달시간 계산
            STT_dgc = GetSTTbrdm(gs, eph, icol, x_dgc(1:3)); % 신호전달시간 계산
            STT_cp = GetSTTbrdm(gs, eph, icol, x_cp(1:3)); % 신호전달시간 계산
            if prn < 300
                if prn > 100 && prn < 200
                    tc = gs - STT ;             % 신호전달시간 보정
                    tc_dgc = gs - STT_dgc ;             % 신호전달시간 보정
                    tc_cp = gs - STT_cp ;             % 신호전달시간 보정
            
                elseif prn > 200 && prn < 300
                    tc = gs - STT- LeapSecBDS;  % 신호전달 시간 보정... bds윤초 반영
                    tc_dgc = gs - STT_dgc- LeapSecBDS;  % 신호전달 시간 보정... bds윤초 반영
                    tc_cp = gs - STT_cp- LeapSecBDS;  % 신호전달 시간 보정... bds윤초 반영
            
                end
                SatPos = GetSatPosNC_GC(eph, icol, tc); % 위성위치 산출
                SatPos = RotSatPos(SatPos, STT); % 지구자전효과 고려
                SatPos_dgc = GetSatPosNC_GC(eph, icol, tc_dgc); % 위성위치 산출
                SatPos_dgc = RotSatPos(SatPos_dgc, STT_dgc); % 지구자전효과 고려
                SatPos_cp = GetSatPosNC_GC(eph, icol, tc_cp); % 위성위치 산출
                SatPos_cp = RotSatPos(SatPos_cp, STT_cp); % 지구자전효과 고려
            
            end
            
            vec_rho = SatPos - vec_site';
            rho = norm(vec_rho);
            [az,el] = xyz2azel(vec_rho, AppLat, AppLon);
            
            Elocal = [sin(az*pi/180)*cos(el*pi/180);...
                cos(az*pi/180)*cos(el*pi/180);...
                sin(el*pi/180)];
            
            Rot = [-sin(AppLon*pi/180) -cos(AppLon*pi/180)*sin(AppLat*pi/180) cos(AppLon*pi/180)*cos(AppLat*pi/180);...
                cos(AppLon*pi/180) -sin(AppLon*pi/180)*sin(AppLat*pi/180) sin(AppLon*pi/180)*cos(AppLat*pi/180);...
                0 cos(AppLat*pi/180) sin(AppLat*pi/180)];
            eECEF = [Rot*Elocal]';
            
            vec_rho_dgc = SatPos_dgc - vec_site_dgc';
            rho_dgc = norm(vec_rho_dgc);
            vec_rho_cp = SatPos_cp - vec_site_cp';
            rho_cp = norm(vec_rho_cp);
            
            
            
            if prn < 200 | (prn > 200 && prn < 300)
                dIono = klobuchar(al, be, gs, az, el, x(1:3)); % 이온층 보정(Klobuchar 모델)
                dTrop = ZHD2SHD(gw, gs, TruePos, el, ZHD); % 대류권 보정
            elseif prn > 300
                dIono = klobuchar(al, be, gs, az, el, x(1:3)); % 이온층 보정(Klobuchar 모델)
                dTrop = ZHD2SHD(gw, gs, TruePos, el, ZHD); % 대류권 보정
            end
            dRel = GetRelBRDC(eph, icol, tc); % 상대성효과
            dRel_dgc = GetRelBRDC(eph, icol, tc_dgc); % 상대성효과
            dtSat = a + b*(tc - toe) + c*(tc - toe)^2 - Tgd + dRel; % 위성시계오차 계산... 그룹딜레이, 상대성효과 보정
            dtSat_dgc = a + b*(tc - toe) + c*(tc - toe)^2 - Tgd + dRel_dgc; % 위성시계오차 계산... 그룹딜레이, 상대성효과 보정
%                         dIono = 0;
%                         dTrop = 0;
            PRC = 0;
            if el >=eleCut % 임계고도각
                W(cnt2,cnt2) = 1;
                %                 W = MakeW_elsnr(el,snr);
                if prn > 100 && prn < 200
                    PRC = GPSPRC(find(GPSPRC(:,1) == gs & GPSPRC(:,2) == prn),3);
                    com = rho + x(4) - CCC*dtSat + dTrop + dIono; % gps 계산값
                    com_dgc = rho_dgc + x_dgc(4) - CCC*dtSat_dgc- PRC; % gps 계산값
                    com_cp = -dTrop - dIono - PRC; % gps 계산값
                    
                    H(cnt2,:) = [-vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 1 0];
                    H_dgc(cnt2,:) = [-vec_rho_dgc(1)/rho_dgc -vec_rho_dgc(2)/rho_dgc -vec_rho_dgc(3)/rho_dgc 1 0];
                    H_cp(cnt2,:) = [-eECEF 1 0];
                    
                    NoGPSsUsed = NoGPSsUsed + 1;
                elseif prn > 200 && prn < 300
                    PRC = BDSPRC(find(BDSPRC(:,1) == gs & BDSPRC(:,2) == prn),3);
                    com = rho + x(5) - CCC*dtSat + dTrop + dIono; % bds 계산값
                    com_dgc = rho_dgc + x_dgc(5) - CCC*dtSat_dgc- PRC; % gps 계산값
                    com_cp = -dTrop - dIono - PRC; % gps 계산값
                    
                    H(cnt2,:) = [-vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 0 1];
                    H_dgc(cnt2,:) = [-vec_rho_dgc(1)/rho_dgc -vec_rho_dgc(2)/rho_dgc -vec_rho_dgc(3)/rho_dgc 0 1];
                    H_cp(cnt2,:) = [-eECEF 0 1];
                    
                    NoBDSsUsed = NoBDSsUsed + 1;
                end
                y(cnt2,:) = obs - com;
                y_dgc(cnt2,:) = obs - com_dgc;
                y_cp(cnt2,:) = com_cp;
                
                cnt2=cnt2+1;
                
            end
        end
        if NoGPSsUsed+NoBDSsUsed >= 5
            HTH = H(1:cnt2-1,:)'*inv(W(1:cnt2-1,1:cnt2-1))*H(1:cnt2-1,:);
            HTy = H(1:cnt2-1,:)'*inv(W(1:cnt2-1,1:cnt2-1))*y(1:cnt2-1,:);
            HTH_dgc = H_dgc(1:cnt2-1,:)'*inv(W(1:cnt2-1,1:cnt2-1))*H_dgc(1:cnt2-1,:);
            HTy_dgc = H_dgc(1:cnt2-1,:)'*inv(W(1:cnt2-1,1:cnt2-1))*y_dgc(1:cnt2-1,:);
            HTH_cp = H_cp(1:cnt2-1,:)'*inv(W(1:cnt2-1,1:cnt2-1))*H_cp(1:cnt2-1,:);
            HTy_cp = H_cp(1:cnt2-1,:)'*inv(W(1:cnt2-1,1:cnt2-1))*y_cp(1:cnt2-1,:);
            
            xhat = inv(HTH) * HTy;
            x = x + xhat;
            
            xhat_dgc = inv(HTH_dgc) * HTy_dgc;
            x_dgc = x_dgc + xhat_dgc;
            xhat_cp = -inv(HTH_cp) * HTy_cp;
            x_cp = x + xhat_cp;                         % DGPS CP alm

            
            if norm(xhat) < EpsStop;
                nEst = nEst + 1;
                estm(nEst,1) = gs;
                estm(nEst,2:6) = x(1:5);
                estm(nEst,7) = NoGPSsUsed;
                estm(nEst,8) = NoBDSsUsed;
                estm(nEst,9) = NoGLOsUsed;
                dgc(nEst,1) = gs;
                dgc(nEst,2:6) = x_dgc(1:5);
                cp(nEst,1) = gs;
                cp(nEst,2:6) = x_cp(1:5);
                
                
%                 fprintf('%8d : %3d : %8.2f : %8.2f : %8.2f\n', j, iter, x(1)' - TruePos(1), x(2)' - TruePos(2), x(3)' - TruePos(3));
                fprintf('%8d : %3d : %8.2f : %8.2f : %8.2f\n', j, iter, x_cp(1)' - TruePos(1), x_cp(2)' - TruePos(2), x_cp(3)' - TruePos(3));
                break;
            end
        else
            break;
        end
    end
end

estm = estm(1:nEst,:);
% Scount=Scount(1:nEst, :);
toc;
[dXYZ, dNEV] = PosTErrorsJOON(estm(:,1), TruePos, estm(:,2:5));
[dXYZ_dgc, dNEV_dgc] = PosTErrorsJOON(dgc(:,1), TruePos, dgc(:,2:5));
[dXYZ_cp, dNEV_cp] = PosTErrorsJOON(cp(:,1), TruePos, cp(:,2:5));
% [dXYZ_cp, dNEV_cp] = PosTErrors4(estm_cp(:,1), TruePos, estm_cp(:,2:5),visiSat);
% [dXYZ_cp_alm, dNEV_cp_alm] = PosTErrors4(estm_cp_alm(:,1), TruePos, estm_cp_alm(:,2:5),visiSat);
[dXYZ_cp_alm, dNEV_cp_alm] = PosTErrorsJOON(cp_alm(:,1), TruePos, cp_alm(:,2:5));
