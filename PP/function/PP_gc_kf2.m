function [estm] = PP_gc_kf(QMfile,eph,TruePos, DOY, YY)

% %
% QMfile = 'QM170125_A';
% QMfile = 'T4_170323_k500';
% QMfile = 'T4_GC';
% QMfile = 'QM170323_Bs';        % 강남 bs 
% QMfile = 'QM170507_Bs';        % 강남1 bs 
% QMfile = 'QM170509_Rv1';        % 강남2 bs 
% QMfile = 'QMfile';
% navfile = 'brdm0250.17p';
% navfile = 'hour0450.17n';
% navfile = 'brdm0820.17p';
% navfile = 'brdm1260.17p';       % 강남 1
% navfile = 'brdm1280.17p';       % 강남 2
% TruePos = [-3041235.578 4053941.677 3859881.013];
% TruePos = [-3052725.48122377,4042747.35562486,3862428.65171442];    % 강남
% TruePos = [-3053365.29481677,4039290.16703344,3865445.80715444];      % 집앞
% YY = 17; DOY =082;
% YY = 17; DOY =126;          % 강남 1
% YY = 17; DOY =128;          % 강남 2
% TruePos = [-3058799.61420451,4083265.35912516,3814946.87192938];

%% 불변 변수 설정: 빛의 속도, 관측치
CCC = 299792458.;   % CCC = Speed of Light [m/s]
%% 임계고도각 설정
eleCut = 15;
%% 윤초 핸들링
LeapSecBDS = 14;

%% 사용할 관측치 선정
g_ObsType = 120; % gps C1
g_ObsType_snr = 141;
g_ObsType_dop = 131;
c_ObsType = 220; % bds C1
c_ObsType_snr = 241;
c_ObsType_dop = 231;

%% QM 파일 읽어들여서 행렬로 저장하고, 사용할 관측치 추출
[arrQM, FinalPRNs, FinalTTs] = ReadQM(QMfile);
arrQM = arrQM(find(arrQM(:,3) < 300),:);
QM = SelectQM_gc(arrQM, g_ObsType, c_ObsType);
QM_snr = SelectQM_gc(arrQM, g_ObsType_snr, c_ObsType_snr);
QM_dop = SelectQM_gc(arrQM, g_ObsType_dop, c_ObsType_dop);
FinalTTs = unique(QM(:,1));
% FinalTTs = FinalTTs(find(FinalTTs(:,1) > 1000),:);

% %% eph 생성
% [eph, trashPRN, trashT]=ReadEPH_all(navfile);

%% 항법메시지를 읽어들여서 행렬로 저장하고, Klobuchar 모델 추출
gps_nav = strcat('brdc', num2str(DOY,'%03d'), '0.', num2str(YY,'%02d'), 'n');   %: Navigation RINEX file
fid = fopen(gps_nav,'r');
if fid == -1
    al = zeros(4,1); be = zeros(4,1);
else
    [al, be] = GetALBE(gps_nav);
end

%% 항법메시지 파일 이름을 이용해 YY, DOY 생성
% YY = str2num(navfile(length(navfile)-2:length(navfile)-1));
% DOY = str2num(navfile(length(navfile)-7:length(navfile)-5));
[gw, gd] = ydoy2gwgd(YY, DOY); %: GPS WEEK 결정

%% 초기좌표 획득
AppPos = TruePos;
gd = xyz2gd(AppPos); AppLat = gd(1); AppLon = gd(2);

%% 추정에 필요한 초기치 설정
Maxiter = 10;
EpsStop = 1e-5;
ctr = 1; ctr2 = 1;
x = [AppPos ctr ctr2]; x = x';

%% 추정과정 시작
NoEpochs = length(FinalTTs);
estm = zeros(NoEpochs, 6);
Scount = zeros(NoEpochs, 1);
% load('PP_170323_gc.mat');
% load('PP_170507_gc.mat');
% load('PP_170509_gc_Bs1.mat');
% load('PP_170509_gc_Rv1.mat');
nEst = 0;

% % P 초기화 
P = eye(5);
P(1:3,1:3) = P(1:3,1:3)*1;
P(4,4) = P(4,4)*0.5; P(5,5) = P(5,5)*0.5; 

% Q 초기화
Q = eye(5);
Q(1,1) = 5;
Q(2,2) = 5;
Q(3,3) = 5;
Q(4,4) = 10000; Q(5,5) = 10000;

for j = 1:NoEpochs-100
    gs = FinalTTs(j);
    indexQM = find(QM(:,1) == gs);
    QM1e = QM(indexQM,:);
    QM1e_snr = QM_snr(indexQM,:);
    QM1e_dop = QM_dop(indexQM,:);
    existprn = intersect(QM1e(:, 2), unique(eph(find(eph(:,22) == 0),18)));
    NoSats = length(existprn);
   
    for iter = 1:1
        HTH = zeros(5,5);
        HTy = zeros(5,1);
        cnt2=1;
        if NoSats <= 6
            break
        end
        vec_site = x(1:3)';
        ZHD = TropGPTh(TruePos, gw, gs); %: TROP: GPT
        for i = 1:NoSats
            NoSatsUsed = 0;
            prn = QM1e(find(QM1e(:,2) == existprn(i)),2);
            obs = QM1e(find(QM1e(:,2) == existprn(i)),4);
            snr = QM1e_snr(find(QM1e_snr(:,2) == existprn(i)),4);
            dop = QM1e_dop(find(QM1e_dop(:,2) == existprn(i)),4);
            icol = PickEPH(eph, prn, gs);
            toe = eph(icol, 8); a = eph(icol, 19); b = eph(icol, 20); c = eph(icol, 21); Tgd = eph(icol, 23);
            %             STT=0;
            STT = GetSTTbrdm(gs, eph, icol, x(1:3)); % 신호전달시간 계산
            if prn > 100 && prn < 200
                tc = gs - STT ;             % 신호전달시간 보정
            elseif prn > 200 && prn < 300
                tc = gs - STT- LeapSecBDS;  % 신호전달 시간 보정... bds윤초 반영
            end
            SatPos = GetSatPosNC_GC(eph, icol, tc); % 위성위치 산출
            SatPos = RotSatPos(SatPos, STT); % 지구자전효과 고려
            vec_rho = SatPos - vec_site';
            rho = norm(vec_rho);
            [az,el] = xyz2azel(vec_rho, AppLat, AppLon);
            if prn > 100 && prn < 200
                %                     dIono = 0;
                %                     dTrop = 0;
                dIono = klobuchar(al, be, gs, az, el, x(1:3)); % 이온층 보정(Klobuchar 모델)
                dTrop = ZHD2SHD(gw, gs, TruePos, el, ZHD); % 대류권 보정
            elseif prn > 200 && prn < 300
                %                     dIono = 0;
                %                     dTrop = 0;
                dIono = klobuchar(al, be, gs, az, el, x(1:3)); % 이온층 보정(Klobuchar 모델)
                dTrop = ZHD2SHD(gw, gs, TruePos, el, ZHD); % 대류권 보정
            end
            %             dRel=0;
            dRel = GetRelBRDC(eph, icol, tc); % 상대성효과
            dtSat = a + b*(tc - toe) + c*(tc - toe)^2 - Tgd + dRel; % 위성시계오차 계산... 그룹딜레이, 상대성효과 보정
            if el >=eleCut % 임계고도각
                W = 1;
%                 W = MakeW_elsnr(el,snr);
                
                elsnr(cnt2,1:6) = [gs, prn, az, el, snr, dop];
                
                if prn > 100 && prn < 200
                    %                         PRC = PickPRC(RTCM,prn,gs);
                    PRC = 0;
                    com = rho + x(4) - CCC*dtSat + dTrop + dIono - PRC; % gps 계산값
                    H(cnt2,:) = [-vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 1 0];
                elseif prn > 200
                    %                         PRC = PickPRC(RTCM,prn,gs);
                    PRC = 0;
                    com = rho + x(5) - CCC*dtSat + dTrop + dIono - PRC; % bds 계산값
                    H(cnt2,:) = [-vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 0 1];
                end
                y(cnt2,1) = obs - com;
                cnt2 = cnt2 + 1;
                NoSatsUsed = NoSatsUsed + 1;
            end
        end
        velocity = GetVelDop(gs, elsnr(:,2), elsnr(:,6), eph, x(1:3)');
        W = PP_cofactor_matrix(elsnr(:,4), elsnr(:,5), 4);
        HTH = H(1:cnt2-1,:)'*inv(W(1:cnt2-1,1:cnt2-1))*H(1:cnt2-1,:);
        HTy = H(1:cnt2-1,:)'*inv(W(1:cnt2-1,1:cnt2-1))*y(1:cnt2-1,:);
        xhat = inv(HTH) * HTy;
        Q(1,1) = abs(velocity(1)); Q(2,2) = abs(velocity(2)); Q(3,3) = abs(velocity(3));
        xp = x;
        Po = P;
        Pp = P + Q;                          %:시스템 노이즈가 공분산에 더해지는 과정
%         K = Pp*H'/(H*Pp*H' + W*1);
        K = Pp*H'/(H*Pp*H' + W*(100/norm(velocity)));
        x = xp + K*y;
        P = Pp - K*H*Pp;
%         P = inv(HTH);
%         xhat = inv(HTH) * HTy;
%         x = x + xhat;
%         if norm(xhat) < EpsStop;
            nEst = nEst + 1;
            estm(nEst,1) = gs;
            estm(nEst,2:6) = x(1:5);
            estm(nEst,7) = NoSatsUsed;
            Scount(nEst,1) = NoSatsUsed;
            VV(nEst,1:5) = [gs, velocity, norm(velocity)];
            fprintf('%8d : %3d : %8.2f : %8.2f\n', j, iter, x(1)' - TruePos(1), x(2)' - TruePos(2));
%             break;
%         end
    end
end
estm = estm(1:nEst,:);
Scount=Scount(1:nEst, :);


for i = 1:length(estm)
    gs = estm(i,1); 
    estm_gd(i,1:4) = [gs xyz2gd(estm(find(estm(:,1) == gs),2:4))];
end

% figure()
% hold on; grid on;
% % plot(Base_gd(:,3), Base_gd(:,2), 'b.')
% % plot(estm_gd(1:2500,3), estm_gd(1:2500,2), 'r.','Markersize',10)
% plot(estm_gd(:,3), estm_gd(:,2), 'r.','Markersize',10)
% plot_google_map