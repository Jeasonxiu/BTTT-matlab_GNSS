%% Double Difference Algorithm using Carrier Phase (Single Frequency)
clear all;close all;%clc;

%% 사용자 설정값 지정
Maximum_Iteration = 1;    % 최대 반복횟수 지정값
stopConditon = 1e-5;      % Iteration 임계값
eleCut = 15;                % 임계고도각 설정

%% 측위에 사용할 QM file 지정.
FileQM_Bs ='QGSDA_17025H1.qm';   %  GSDA (가산동 옥상 A)
FileQM_Rv ='QGSDB_17025H1.qm';   %  GSDB (가산동 옥상 B)
%% Part 1. Between-Satellite Single Difference Positioning algorithms.

CCC = 299792458.;           % 빛의 속도 상수 선언. [m/s]
Freq_L1 = 1575.42e6;        % L1 주파수
Lambda_L1 = CCC/Freq_L1;    % L1 파장

[arrQM_Bs, FinalPRNs_Bs, FinalTTs_Bs] = ReadQM(FileQM_Bs);  % 측위에 사용할 QM file을 읽어들임.
[arrQM_Rv, FinalPRNs_Rv, FinalTTs_Rv] = ReadQM(FileQM_Rv);  % 측위에 사용할 QM file을 읽어들임.
FinalTTs = unique(intersect(FinalTTs_Bs,FinalTTs_Rv));      % 두 수신 데이터의 공통 수신시각 추출

ObsType = 111;                                  % 관측치 타입 :  L1 (CP)
ObsType_CA = 120; % 시간지연 보정으로 코드 관측치 사용 - 120: CA

ObsType_Selected_QM_Bs = SelectQM(arrQM_Bs, ObsType);       % 측위에 사용할 관측치만 추출. (in QM_Bs)
ObsType_Selected_QM_Rv = SelectQM(arrQM_Rv, ObsType);       % 측위에 사용할 관측치만 추출. (in QM_Rv)

% 해당 QM file에 상응하는 brdc 데이터명 지정.
YY = FileQM_Bs(7:8);
DOY = FileQM_Bs(9:11);
FileNav = strcat('brdc',DOY,'0.',YY,'n');

eph = ReadEPH(FileNav); % 해당 brdc 데이터를 읽어들임.

[gw, gd] = ydoy2gwgd(str2num(YY), str2num(DOY)); % 계산에 이용할 GPS week number 산출.

% 사용자 좌표
User_name_Bs = FileQM_Bs(2:5);
User_name_Rv = FileQM_Rv(2:5);
switch (User_name_Bs)
    case 'GSDA'
        TruePos_Bs = [-3041235.578 4053941.677 3859881.013];       % GSDA   
    case 'GSDB'
        TruePos_Bs = [-3041241.741 4053944.143 3859873.640];       % GSDB            
end

switch (User_name_Rv)
    case 'GSDA'
        TruePos_Rv = [-3041235.578 4053941.677 3859881.013];       % GSDA   
    case 'GSDB'
        TruePos_Rv = [-3041241.741 4053944.143 3859873.640];       % GSDB          
end
AppPos_Rv = TruePos_Rv+10;

%% Kalman Filter를 이용한 측위 과정
% 추정을 위한 변수값 설정

NoEpochs = length(FinalTTs); % Number of Epochs.

Estimated_Position = zeros(NoEpochs, 5);   % 추정치 저장을 위한 행렬. (열 순서대로 시간, X, Y, Z, 시계오차)
No_Est = 0;

tic
%% 칼만설정
    Sats_for_Process = FinalPRNs_Bs;
    NoSats = length(Sats_for_Process);
    
    %%
    x = zeros(3+(NoSats-1),1); % 측위 결과 계산을 위한 초기값 지정.
    x(1:3) = AppPos_Rv';
        
    A = eye(3+(NoSats-1));
    A = 1 .* A;

    Q = eye(3+(NoSats-1));
    P = eye(3+(NoSats-1));

    save_P = zeros(1,(length(P))+2);
%% 고고씽
for i = 1:NoEpochs
    %%
    IndexS_QM_Bs = find(ObsType_Selected_QM_Bs(:,1) == FinalTTs(i));   % 각 시각의 L1 measurements 만을 모으기 위한 index.
    IndexS_QM_Rv = find(ObsType_Selected_QM_Rv(:,1) == FinalTTs(i));   % 각 시각의 L1 measurements 만을 모으기 위한 index.
    
    sub_QM_Bs_temp2 = ObsType_Selected_QM_Bs(IndexS_QM_Bs,:);
    sub_QM_Rv_temp2 = ObsType_Selected_QM_Rv(IndexS_QM_Rv,:);
    
    Index_QM_Bs = NaN(1,1);
    Index_QM_Rv = NaN(1,1);
    
    for abcd = 1:NoSats
        Index_QM_Bs_temp2 = find(sub_QM_Bs_temp2(:,2) == Sats_for_Process(abcd));
        Index_QM_Bs = union(Index_QM_Bs,Index_QM_Bs_temp2);
        Index_QM_Rv_temp2 = find(sub_QM_Rv_temp2(:,2) == Sats_for_Process(abcd));
        Index_QM_Rv = union(Index_QM_Rv,Index_QM_Rv_temp2);
    end
    
    sub_QM_Bs = sub_QM_Bs_temp2(Index_QM_Bs(1:NoSats),:);
    sub_QM_Rv = sub_QM_Rv_temp2(Index_QM_Rv(1:NoSats),:);
    
    gs = FinalTTs(i);               % 계산에 이용할 GPS second 추출.
    
    rcvPos = AppPos_Rv;                % 수신기 xyz 좌표
    
    %% 기준위성(RS) 선별 //  SatsAZEI - c1(gs), c2(prn), c3(az), c4(el)
    [SatsAZEl, indxRS]=pickRSel(gs,Sats_for_Process, eph, rcvPos, Correction);      % RS : Reference Satellite 기준위성
    RS_PRN = Sats_for_Process(indxRS);
    
    %% 기준위성(RS)의 관측치 추출
    obs_RS_Bs = sub_QM_Bs(find(sub_QM_Bs(:,2)==RS_PRN),4);
    obs_RS_Rv = sub_QM_Rv(find(sub_QM_Rv(:,2)==RS_PRN),4);
    
    %% 반복 계산 시작
        
        % 최소자승법 계산을 위한 Matrix 선언.
        H = zeros(1,3+(NoSats-1));
        y = zeros(1,1);             % y(관측치-계산값) 선언
        
        R = eye(NoSats-1,NoSats-1);  % R matrix 선언
        
        % 반복계산때마다 새로운 x값을 사용
        vec_site = x(1:3)';              % 수신기 xyz 좌표 (반복계산이 끝나기 전까지 갱신됨.)
        
        %% 기준위성 좌표 먼저 계산
        % 해당 시각에 근접한 데이터를 이용하기 위한 column 위치 산출.
        indxEPH_RS = PickEPH(eph, RS_PRN, gs);
        if indxEPH_RS == 0
            continue;
        end
        
        STT = obs_RS_Bs(i)/CCC;

        tc = gs - STT;
        vec_Sat_RS = GetSatPos_G(eph, indxEPH_RS, tc);      % 해당 위성의 위치 산출
        
        vec_Sat_RS = RotSatPos(vec_Sat_RS, STT);
        
        vec_rho_RS_Bs = vec_Sat_RS - TruePos_Bs; 
        vec_rho_RS_Rv = vec_Sat_RS - vec_site; 
        rho_RS_Bs = norm(vec_rho_RS_Bs);
        rho_RS_Rv = norm(vec_rho_RS_Rv);

        No_OS = 0;
        
        %% 각 위성에 대한 관측치, 계산치, H행렬계산
        for j = 1:NoSats

            el_OS = SatsAZEl(j,4);

            OS_PRN = Sats_for_Process(j);         	  % 해당 위성의 PRN
            
            if OS_PRN == RS_PRN
                continue
            else
                No_OS = No_OS+1;
            end
            
            %% DD관측치 생성파트
            obs_OS_Bs = sub_QM_Bs(find((sub_QM_Bs(:,1)==gs)&(sub_QM_Bs(:,2)==OS_PRN)),4);    % y 행렬 계산을 위한 해당 시각, 해당 위성의 관측값(obs) 추출.
            obs_OS_Rv = sub_QM_Rv(find((sub_QM_Rv(:,1)==gs)&(sub_QM_Rv(:,2)==OS_PRN)),4);    % y 행렬 계산을 위한 해당 시각, 해당 위성의 관측값(obs) 추출.

            % 해당 시각에 근접한 데이터를 이용하기 위한 column 위치 산출.
            indxEPH_OS = PickEPH(eph, OS_PRN, gs);   
            if indxEPH_OS == 0
                continue;
            end

            STT = obs_OS_Bs/CCC;

            tc = gs - STT;
            vec_Sat_OS = GetSatPos_G(eph, indxEPH_OS, tc);       % 해당 위성의 위치 산출
            
            vec_Sat_OS = RotSatPos(vec_Sat_OS, STT);
        
            %% SD 계산치 생성파트-거리 계산치를 각각 계산한 다음 DD계산치 계산

            vec_rho_OS_Bs = vec_Sat_OS - TruePos_Bs;     % 위성과 수신기간의 rho vector 산출.
            vec_rho_OS_Rv = vec_Sat_OS - vec_site;     % 위성과 수신기간의 rho vector 산출.
            rho_OS_Bs = norm(vec_rho_OS_Bs);              % 계산 결과에 따른 위성으로부터 수신기까지의 거리 계산.
            rho_OS_Rv = norm(vec_rho_OS_Rv);              % 계산 결과에 따른 위성으로부터 수신기까지의 거리 계산.

            com = ((rho_RS_Bs - rho_OS_Bs) - (rho_RS_Rv - rho_OS_Rv));
            obs = ((obs_RS_Bs - obs_OS_Bs) - (obs_RS_Rv - obs_OS_Rv));
            
            if el_OS > eleCut %% - 각도를 바꿔가면서 해보고
                R(No_OS,No_OS) = 1;
            else
                R(No_OS,No_OS) = ( 1/sind(el_OS) );
            end
                
            y(No_OS,1) = obs - com;    % 관측치 - 계산값
            
            % H Matrix 산출.
            for k = 1:3
                H(No_OS,k) = ( (vec_rho_RS_Rv(k)) / rho_RS_Rv ) - ( (vec_rho_OS_Rv(k)) / rho_OS_Rv );
            end
            
        end

        xp = A*x;
        Pp = A*P*A' + Q;
        K = Pp*H'*inv(H*Pp*H' + R);

        x = xp + K*(y);
        P = Pp - K*H*Pp;

        No_Est = No_Est + 1;
        Estimated_Position(No_Est,1) = gs;            % 추정치 저장. (행 순서대로 시간, X, Y, Z, 시계오차)
        Estimated_Position(No_Est,2:5) = x(1:4);      % 추정치 저장. (행 순서대로 시간, X, Y, Z, 시계오차)
 
        for ppppp = 1:(length(P))
            save_P(No_Est,ppppp) = P(ppppp,ppppp);
        end
        save_P(No_Est,((length(P))+1)) = sub_QM_Bs(1,1);
        save_P(No_Est,((length(P))+2)) = NoSats;

end
toc


%% Part 2. Positioning Accuracy Plotting 
Estimated_Position = Estimated_Position(1:No_Est, :);
estmTT = Estimated_Position(:, 1);        % 측위 결과에 상응하는 Time Tags (gps second)
estmPos = Estimated_Position(:, 2:4);     % 추정 계산된 수신기위치
NoPos = length(estmPos);    % 추정 결과 데이터의 길이.

geodetic_true = xyz2gd(TruePos_Rv); % TRUE POSITION Latitude, Longitude
TrueLat = geodetic_true(1);      % TRUE POSITION Latitude
TrueLon = geodetic_true(2);      % TRUE POSITION Longitude

% delta X, delta Y, delta Z 산출
dXYZ = zeros(NoPos,3);
for k = 1:NoPos
    dXYZ(k,:) = estmPos(k,:) - TruePos_Rv;
end

dNEV = xyz2topo(dXYZ, TrueLat, TrueLon); % delta XYZ를 delta NEV로 변환


% -------------------------------------------------------------------------
% RMSE 산출
dN = dNEV(:,1); 
dE = dNEV(:,2); 
dNE = sqrt(dN.^2 + dE.^2);       
rmsH = sqrt(mean(dNE.^2));  % 수평 오차 산출
Last_H = sqrt(dNE(end)^2);  % 수평 오차 산출

mean_dN = sqrt(mean(dN.^2));            % 남북방향 수평 오차 산출
Last_dN = sqrt(dN(end)^2);            % 남북방향 수평 오차 산출
mean_dE = sqrt(mean(dE.^2));            % 동서방향 수평 오차 산출
Last_dE = sqrt(dE(end)^2);            % 동서방향 수평 오차 산출

dV = dNEV(:,3);
rmsV = sqrt(mean(dV.^2));   % 수직 오차 산출
Last_V = sqrt(dV(end)^2);   % 수직 오차 산출

d3 = sqrt(dN.^2 + dE.^2 + dV.^2); 
rms3 = sqrt(mean(d3.^2));   % 3차원 오차 산출
Last_3 = sqrt(d3(end)^2);   % 3차원 오차 산출


fprintf('\nRMSE Values \n  H :%8.2f      ( N :%5.2f / E :%5.2f )\n  V :%8.2f \n 3-D:%8.2f\n', rmsH, mean_dN, mean_dE, rmsV, rms3)
fprintf('\nRMSE Values (Last Point) \n  H :%8.2f      ( N :%5.2f / E :%5.2f )\n  V :%8.2f \n 3-D:%8.2f\n', Last_H, Last_dN, Last_dE, Last_V, Last_3)

% -------------------------------------------------------------------------
    % Plotting (delta H)
    rXY = max([abs(min(dNEV(:,1))) abs(max(dNEV(:,1))) abs(min(dNEV(:,2))) abs(max(dNEV(:,2)))]);
    rXY = ceil(rXY);
    
    positionVector1 = [0.04,0.36,0.34,0.34];
    subplot('Position',positionVector1)
    set(gca,'FontSize',13,'FontWeight','Bold')
    hold on;
    plot(dE, dN,'o');
    plot(dE(end), dN(end),'ro','linewidth',3,'markersize',10);

    axis([-rXY rXY -rXY rXY]); grid on;
    xlabel('\Delta E (meters)'); ylabel('\Delta N (meters)')
    
    % Plotting (delta N, delta E)
    tHour = mod(estmTT, 86400);
    tHour = tHour/3600;          % 출력 시간단위를 Hour로 설정.
    
    positionVector2 = [0.42,0.70,0.56,0.28];
    subplot('Position',positionVector2)
    hold on;
    plot(tHour, dN, '.r:', tHour, dE, '.b:'); axis([min(tHour) max(tHour) min(min(dN),min(dE)) max(max(dN),max(dE))]); grid on;
    legend('\Delta N', '\Delta E')
    set(gca,'FontSize',13,'FontWeight','Bold')
    ylabel('\Delta H (meters)');
    
    % Plotting (delta V)
    positionVector3 = [0.42,0.38,0.56,0.28];
    subplot('Position',positionVector3)
    hold on;
    plot(tHour, dV, '.:'); axis([min(tHour) max(tHour) min(0,min(dV)) max(dV)]); grid on;
    set(gca,'FontSize',13,'FontWeight','Bold')
    ylabel('\Delta U (meters)')
    
    % Plotting (delta 3-D)
    positionVector4 = [0.42,0.06,0.56,0.28];
    subplot('Position',positionVector4)
    hold on;
    plot(tHour, d3, '.:'); axis([min(tHour) max(tHour) min(0,min(d3)) max(d3)]); grid on;
    set(gca,'FontSize',13,'FontWeight','Bold')
    ylabel('\Delta 3-D (meters)')
    xlabel('Hours');
   
    figure(2)
    save_P(:,((length(P))+1)) = mod(save_P(:,((length(P))+1)),86400)/3600;
    for plotp = 1:length(P)
        hold on;
        
        if (35+plotp) > 42
            NoStyle = mod((35+plotp),42) + 35;
        else
            NoStyle = 35+plotp;
        end
            
        Plot_XY(save_P(:,((length(P))+1)),save_P(:,plotp),NoStyle,10,10)
    end

    set(gca,'FontSize',20,'FontWeight','Bold')
    ylabel('Elements of Matrix P')
    xlabel('Hours');    
    grid on;
       
    
