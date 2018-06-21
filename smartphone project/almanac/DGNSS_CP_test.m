clear all
close all

% %% nmea로깅 파일
% % filename = 'S7170912.txt'; year = 2017; month = 09; day = 12;   % S7 
% filename = 'N9170914.txt'; year = 2017; month = 09; day = 14;   % NEXUS9 
% DOY = date2doy(day,month,year);
% %% nmea 행렬 생성
% [nmea] = writeNMEA(filename,year,month,day);
% %% 로깅시간
% FinalTTs = unique(nmea(:,1));
% 
% % TruePos = [-3041235.578 4053941.677 3859881.013];   % : JPspace A point
% TruePos = [-3041241.741 4053944.143 3859873.640];   % : JPspace B point
% 
% %% load a PRC file
% % PRCfile = 'PPS1_170912.t41';
% PRCfile = 'PPS1_170914.t41';
% %% 로깅시간 내 PRC정보 sorting
% [GPSPRC, BDSPRC, GLOPRC, PRC_sorted] = PRCsort(PRCfile, FinalTTs);
% %% ZHD 생성을 위한 GPS week 생성
% [gw, GD] = ydoy2gwgd(year-2000, DOY); %: GPS WEEK 결정
% %% 항법메시지를 읽어들여서 행렬로 저장하고, Klobuchar 모델 추출
% gps_nav = strcat('brdc', num2str(DOY,'%03d'), '0.', num2str(year-2000,'%02d'), 'n');   %: Navigation RINEX file
% fid = fopen(gps_nav,'r');
% if fid == -1
%     al = zeros(4,1); be = zeros(4,1);
% else
%     [al, be] = GetALBE(gps_nav);
% end


% load('DGNSS_CP_170912.mat');
load('DGNSS_CP_171016_2.mat');
% load('DGNSS_CP_170914.mat');
TruePos = [-3041235.578 4053941.677 3859881.013];   % : JPspace A point
FinalTTs = intersect(FinalTTs, PRC_sorted(:,1));
% PRC_sorted = PRC_sorted(find(PRC_sorted(:,2) < 200),:);
% nmea = nmea(find(nmea(:,7) < 200),:);
for i=1:length(FinalTTs)
    %% 1 Epoch 행렬 생성
    gs = FinalTTs(i);
    RcvLLH(i,:) = nmea(min(find(nmea(:,1) == gs)),1:4);
    RcvXYZ(i,:) = [gs gd2xyz(RcvLLH(i,2:4))];
    Epoch = nmea(find(nmea(:,1) == gs),:);
    Epoch = Epoch(find(Epoch(:,8) ~= 0 & Epoch(:,9) ~= 0),:);
    %% PRC 1Epoch
    PRC1e = PRC_sorted(find(PRC_sorted(:,1) == gs),:);
    ExistPRN = intersect(Epoch(:,7), PRC1e(:,2));
    %% 현재 수신된 Epoch의 수가 계산하기에 부족할 경우
    if length(ExistPRN) < 5
        break;
    end
    %% 사용 시스템 확인
    NumGPS = length(Epoch(find(Epoch(:,7) < 200)));
    if isempty(NumGPS)
        NumGPS = 0; end
    NumGLO = length(Epoch(find(Epoch(:,7) > 300)));
    if isempty(NumGLO)
        NumGLO = 0; end
    NumBDS = length(Epoch(find(Epoch(:,7) > 200 & Epoch(:,7) < 300)));
    if isempty(NumBDS)
        NumBDS = 0; end
    if NumGPS ~= 0 & NumGLO ~= 0 & NumBDS ~= 0
        EstPos = [RcvXYZ(i,2:4), 1, 1, 1]';
    elseif NumGPS ~= 0 & NumGLO ~= 0
        EstPos = [RcvXYZ(i,2:4), 1, 1]';
    elseif NumGPS ~= 0 & NumBDS ~= 0
        EstPos = [RcvXYZ(i,2:4), 1, 1]';
    elseif NumGLO ~= 0 & NumBDS ~= 0
        EstPos = [RcvXYZ(i,2:4), 1, 1]';
    elseif NumGPS ~= 0
        EstPos = [RcvXYZ(i,2:4), 1]';
    elseif NumGLO ~= 0
        EstPos = [RcvXYZ(i,2:4), 1]';
    elseif NumBDS ~= 0
        EstPos = [RcvXYZ(i,2:4), 1]';
    end
    %% ZHD 계산
    ZHD = TropGPTh(EstPos, gw, gs); %: TROP: GPT
    %% Iteration 시작
    Maxiter = 1;
    for iter = 1:Maxiter
        %% 관측행렬 생성을 위한 카운트
        cnt2=1;
        for j=1:length(ExistPRN)
            %% NMEA 해당 위성 선택
            prn = ExistPRN(j);
            %% 해당위성의 az, el
            az = Epoch(find(Epoch(:,7) == prn),9); el = Epoch(find(Epoch(:,7) == prn),8);
            %% NMEA az, el 정보로부터 Local frame 생성
            Elocal = [sin(az*pi/180)*cos(el*pi/180);...
                cos(az*pi/180)*cos(el*pi/180);...
                sin(el*pi/180)];
            %% 현재 gs에서의 수신기 Lat, Lon
            AppLLH = xyz2gd([EstPos(1:3)]');
            AppLat = AppLLH(1); AppLon = AppLLH(2);
            %% 현재 gs에서의 회전행렬 생성
            Rot = [-sin(AppLon*pi/180) -cos(AppLon*pi/180)*sin(AppLat*pi/180) cos(AppLon*pi/180)*cos(AppLat*pi/180);...
                cos(AppLon*pi/180) -sin(AppLon*pi/180)*sin(AppLat*pi/180) sin(AppLon*pi/180)*cos(AppLat*pi/180);...
                0 cos(AppLat*pi/180) sin(AppLat*pi/180)];
            %% 시선백터 생성
            eECEF = [Rot*Elocal]';
            dIono = klobuchar(al, be, gs, az, el, EstPos(1:3)); % 이온층 보정(Klobuchar 모델)
            dTrop = ZHD2SHD(gw, gs, EstPos(1:3), el, ZHD); % 대류권 보정
                        
            if NumGPS ~= 0 & NumGLO ~= 0 & NumBDS ~= 0
                
                if prn < 200
                    %% 해당위성의 PRC
                    PRC = GPSPRC(find(GPSPRC(:,1) == gs & GPSPRC(:,2) == prn),3);
                    H(cnt2,:) = [-eECEF, 1, 0, 0];
                elseif prn > 300
                    PRC = GLOPRC(find(GLOPRC(:,1) == gs & GLOPRC(:,2) == prn),3);
                    H(cnt2,:) = [-eECEF, 0, 1, 0];
                else
                    PRC = BDSPRC(find(BDSPRC(:,1) == gs & BDSPRC(:,2) == prn),3);
                    H(cnt2,:) = [-eECEF, 0, 0, 1];
                end
            elseif NumGPS ~= 0 & NumGLO ~= 0
                if prn < 200
                    %% 해당위성의 PRC
                    PRC = GPSPRC(find(GPSPRC(:,1) == gs & GPSPRC(:,2) == prn),3);
                    H(cnt2,:) = [-eECEF, 1, 0];
                elseif prn > 300
                    PRC = GLOPRC(find(GLOPRC(:,1) == gs & GLOPRC(:,2) == prn),3);
                    H(cnt2,:) = [-eECEF, 0, 1];
                end
            elseif NumGPS ~= 0 & NumBDS ~= 0
                if prn < 200
                    %% 해당위성의 PRC
                    PRC = GPSPRC(find(GPSPRC(:,1) == gs & GPSPRC(:,2) == prn),3);
                    H(cnt2,:) = [-eECEF, 1, 0];
                elseif prn > 200 & prn < 300
                    PRC = BDSPRC(find(BDSPRC(:,1) == gs & BDSPRC(:,2) == prn),3);
                    H(cnt2,:) = [-eECEF, 0, 1];
                end
            elseif NumGLO ~= 0 & NumBDS ~= 0
                if prn > 300
                    %% 해당위성의 PRC
                    PRC = GLOPRC(find(GLOPRC(:,1) == gs & GLOPRC(:,2) == prn),3);
                    H(cnt2,:) = [-eECEF, 1, 0];
                elseif prn > 200 & prn < 300
                    PRC = BDSPRC(find(BDSPRC(:,1) == gs & BDSPRC(:,2) == prn),3);
                    H(cnt2,:) = [-eECEF, 0, 1];
                end
            elseif NumGPS ~= 0
                PRC = GPSPRC(find(GPSPRC(:,1) == gs & GPSPRC(:,2) == prn),3);
                H(cnt2,:) = [-eECEF, 1];
            elseif NumGLO ~= 0
                PRC = GLOPRC(find(GLOPRC(:,1) == gs & GLOPRC(:,2) == prn),3);
                H(cnt2,:) = [-eECEF, 1];
            elseif NumBDS ~= 0
                PRC = BDSPRC(find(BDSPRC(:,1) == gs & BDSPRC(:,2) == prn),3);
                H(cnt2,:) = [-eECEF, 1];
            end
            W(cnt2,cnt2) = 1;
%             Y(cnt2,:) = -dIono -dTrop - PRC;
%             Y(cnt2,:) = -dIono - dTrop;   % best1
            Y(cnt2,:) = -dIono;   % S7171016
%             Y(cnt2,:) = PRC;   % S7170914
%             Y(cnt2,:) = +dIono*1 + dTrop*1 +PRC;   % 
            cnt2=cnt2+1;
        end
        HTH = H(1:cnt2-1,:)'*inv(W(1:cnt2-1,1:cnt2-1))*H(1:cnt2-1,:);
        HTy = H(1:cnt2-1,:)'*inv(W(1:cnt2-1,1:cnt2-1))*Y(1:cnt2-1,:);
        xhat = -inv(HTH) * HTy;
        XHAT(i,:) = [gs xhat(1:3)'];
        EstPos = EstPos + xhat;
%         norm(xhat)
%         if norm(xhat) < 1e-5
            if iter == Maxiter
            estm(i,1) = gs;
            estm(i,2:4) = EstPos(1:3);
            fprintf('%8d : %3d : %8.2f : %8.2f : %8.2f\n', i, iter, EstPos(1)' - TruePos(1), EstPos(2)' - TruePos(2), EstPos(3)' - TruePos(3));
%             break;
        end
    end
    EstPos = [RcvXYZ(i,2:4), 1, 1, 1]';
end
[dXYZ, dNEV] = PosTErrorsJOON(estm(:,1), TruePos, estm(:,2:4));
[dXYZ, dNEV] = PosTErrorsJOON(RcvXYZ(:,1), TruePos, RcvXYZ(:,2:4));
