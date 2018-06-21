clear all;
close all;

%% Base, Rover 관측데이터
% 관측 YY, DOY
% YY = '18'; DOY = '073';     % JAVAD-JAVAD
YY = '18'; DOY = '107';     % JAVAD-JAVAD
% 관측데이터 로딩
% QMBs = load('QDSAJ_18073m');    % 18073 JAVAD-A
% QMRv = load('QDSBJ_18073m');    % 18073 JAVAD-B
QMBs = load('QDSAJ_18107m');    % 18107 JAVAD-A
QMRv = load('QDSBJ_18107m');    % 18107 JAVAD-Bsscanf

QMBs_C1 = QMBs(QMBs(:,3) == 120, :);        % Base C1 데이터 추출
QMBs_L1 = QMBs(QMBs(:,3) == 111, :);        % Base L1 데이터 추출
QMRv_C1 = QMRv(QMRv(:,3) == 120, :);        % Base C1 데이터 추출
QMRv_L1 = QMRv(QMRv(:,3) == 111, :);        % Base L1 데이터 추출
% 관측데이터 기반 기본 파라미터 선언
% NoSats = 7; prn_RS = 16;        % 현재 초기버전에서는 고정 18073 javad
NoSats = 6; prn_RS = 14;        % 현재 초기버전에서는 고정 18073 javad
FinalTTs = intersect(unique(QMBs_L1(:,1)), unique(QMRv_L1(:,1)));
FinalTT = FinalTTs(1:3988); FinalTTs = FinalTT;
%% 항법메시지 로딩
FileNav = strcat('brdc',DOY,'0.',YY,'n');   % GPS 방송궤도력 파일 읽어들이기 =>  'brdc0000.00n'
eph = ReadEPH(FileNav); % 해당 brdc 데이터를 읽어들임.

%% 포지셔닝 기본 파라미터 선언
CCC = 299792458.;           	% 빛의 속도 상수 선언. [m/s]
Freq_L1 = 1575.42e6;    Lambda = CCC/Freq_L1;    % GPS L1신호 주파수/파장

%% Base, Rover TruePos 선언
TruePos_Bs = [-3041235.578 4053941.677 3859881.013];        % A Point
TruePos = [-3041241.741 4053944.143 3859873.640];           % B point
TruePos_gd = xyz2gd(TruePos);                               % TruePos geodetic coordinates
%% Kalman-filter 초기값 선언
Q = MakeQkt(NoSats-1);
P = eye(3+NoSats-1);                            % xyz + 기준위성을 제외한 NoSats
P(1:3,1:3) = P(1:3,1:3) * 2;                    % P xyz 초기값 선언
P(4:3+NoSats-1,4:3+NoSats-1) = P(4:3+NoSats-1,4:3+NoSats-1) * 10;       % P 위성 별 N 초기값 선언
A = eye(3+NoSats-1);                            % xyz + 기준위성을 제외한 NoSats
R = eye(NoSats-1,NoSats-1);                     % 기준위성과 나머지 위성 조합수
%% Rover 초기좌표
AppPos = TruePos + 1;               % 기준좌표에 1m 오차 더함
% 이중차분 연산 시 로버 좌표 갱신을 위한 x 값 선언
x = AppPos';

%% GPS-RTK 이중차분 시작
% 이중차분 초기값 선언
N = zeros(NoSats-1,1);
H = zeros(1,3+NoSats-1);
y = zeros(NoSats-1,1);

for i = 1:length(FinalTTs)
    vec_site = x(1:3)';
    gs = FinalTTs(i);                               % 현재 Epoch gs
    % Base, Rover 현재 Epoch 관측데이터 추출
    QMBs_C1_1e = QMBs_C1(QMBs_C1(:,1) == gs, :);    % 현재 gs의 C1 관측값
    QMBs_L1_1e = QMBs_L1(QMBs_L1(:,1) == gs, :);    % 현재 gs의 L1 관측값
    QMRv_C1_1e = QMRv_C1(QMRv_C1(:,1) == gs, :);    % 현재 gs의 C1 관측값
    QMRv_L1_1e = QMRv_L1(QMRv_L1(:,1) == gs, :);    % 현재 gs의 L1 관측값
    SatsList = intersect(QMBs_L1_1e(:,2),QMRv_L1_1e(:,2));
    
    % 이중차분을 위한 관측치, 위성 좌표 및, 방위각, 고도각 산출
    for j = 1:NoSats
        prn = SatsList(j);
        % 현재 gs와 prn에 맞는 eph column 선택
        icol = PickEPH(eph,prn,gs);
        toe = eph(icol, 8); a = eph(icol, 19); b = eph(icol, 20); c = eph(icol, 21); Tgd = eph(icol, 23);
        %% Base 관측정보 및 위성 정보 행렬 생성
        % Base 관측정보 행렬 생성
        C1_Bs = QMBs_C1_1e(QMBs_C1_1e(:,2) == prn,4);      % C1 observation value(Base)
        L1_Bs = QMBs_L1_1e(QMBs_L1_1e(:,2) == prn,4);      % L1 observation value(Base)
        Obs_Bs(j,:) = [prn, C1_Bs, L1_Bs];
        % Base 위성 위치 계산
        STT_Bs = GetSTTbrdc(gs, prn, eph, TruePos_Bs);
        tc_Bs = gs - STT_Bs;
        vec_sat_Bs = GetSatPosNC(eph, icol, tc_Bs);
        vec_sat_Bs = RotSatPos(vec_sat_Bs, STT_Bs);      % 위성 좌표
        vec_rho_Bs = vec_sat_Bs - TruePos_Bs';
        rho_Bs = norm(vec_rho_Bs);                    % Rho 연산
        % Base 위치에서의 위성 방위각, 고도각 계산
        Bs_gd = xyz2gd(TruePos_Bs);             % 위성의 고도각 계산을 위한 Base 좌표 geodetic 변환
        [az,el] = xyz2azel(vec_rho_Bs,Bs_gd(1),Bs_gd(2));
        SatsInfo_Bs(j,:) = [prn, vec_rho_Bs', rho_Bs, az, el];
        %% Rover 관측정보 및 위성 정보 행렬 생성
        % Rover 관측정보 행렬 생성
        C1_Rv = QMRv_C1_1e(QMRv_C1_1e(:,2) == prn,4);      % C1 observation value(Rover)
        L1_Rv = QMRv_L1_1e(QMRv_L1_1e(:,2) == prn,4);      % L1 observation value(Rover)
        Obs_Rv(j,:) = [prn, C1_Rv, L1_Rv];
        % Rover 위성 위치 계산
        STT = GetSTTbrdc(gs, prn, eph, vec_site);
        tc = gs - STT;
        vec_sat = GetSatPosNC(eph, icol, tc);
        vec_sat = RotSatPos(vec_sat, STT);      % 위성 좌표
        vec_rho = vec_sat - vec_site';
        rho = norm(vec_rho);                    % Rho 연산
        % Rover 현재위치에서의 위성 방위각, 고도각 계산
        Rv_gd = xyz2gd(vec_site);             % 위성의 고도각 계산을 위한 Rover 좌표 geodetic 변환
        [az,el] = xyz2azel(vec_rho,Rv_gd(1),Rv_gd(2));
        SatsInfo_Rv(j,:) = [prn, vec_rho', rho, az, el];
    end
    
    %% 기준위성 선정
    % 현재는 정해놓은 위성을 사용
    idx_RS = find(SatsList(:) == prn_RS);
    SatsList_OS = SatsList(SatsList(:) ~= prn_RS);         % 기준위성을 제외한 위성목록
    
    %% 모호정수 초기값 생성
    if N(1) == 0
        N = N_init(Obs_Bs, Obs_Rv, idx_RS);             % 첫 Epoch 에서 N 초기값 생성
    else
        N = x(4:end);
    end
    %% R: measurement noise
    R = MakeRkt_(SatsInfo_Rv);     % 위성 고도각을 이용한 R 어레이 생성 : 김준성
%     R = MakeRkt(SatsInfo_Rv);     % 위성 고도각을 이용한 R 어레이 생성 : 교수님
    % 기준위성 제외
    R(:,idx_RS) = []; R(idx_RS,:) = [];
    % H matrix 초기화
    H=[];
    %% 이중차분 계산값 산출
    for k=1:length(SatsList_OS)
        prn_OS = SatsList_OS(k);
        idx_OS = find(SatsInfo_Rv(:,1) == prn_OS);
        x(3+k) = N(k);
        if idx_OS == idx_RS
            continue;
        end
        % com 산출
        com = (SatsInfo_Bs(idx_RS,5)-SatsInfo_Rv(idx_RS,5)) -...
            (SatsInfo_Bs(idx_OS,5)-SatsInfo_Rv(idx_OS,5)) +...
            Lambda*(N(k));
        COM(k,1) =com;
        obs = ((Obs_Bs(idx_RS,3)-Obs_Rv(idx_RS,3)) - (Obs_Bs(idx_OS,3)-Obs_Rv(idx_OS,3)))*Lambda;
        y(k,1) = obs -com;
        % H 행렬 생성
        H(k,1) = (SatsInfo_Rv(idx_RS,2)/SatsInfo_Rv(idx_RS,5)) - (SatsInfo_Rv(idx_OS,2)/SatsInfo_Rv(idx_OS,5));
        H(k,2) = (SatsInfo_Rv(idx_RS,3)/SatsInfo_Rv(idx_RS,5)) - (SatsInfo_Rv(idx_OS,3)/SatsInfo_Rv(idx_OS,5));
        H(k,3) = (SatsInfo_Rv(idx_RS,4)/SatsInfo_Rv(idx_RS,5)) - (SatsInfo_Rv(idx_OS,4)/SatsInfo_Rv(idx_OS,5));
        % 모호정수 편미분값
        H(k, 3+k) = Lambda;
    end
%     H(idx_RS,:) =[]; y(idx_RS) =[];
    %% Kalman-filter 연산
    xp = A*x;
    Pp = A*P*A' + Q;
    K = Pp*H'/(H*Pp*H' + R);
    % 좌표 갱신
    x = xp + K*(y);
    % 공분산 행렬 갱신
    P = Pp - K*H*Pp;
    estm(i,:) = [gs, x'];
    %% dXYZ, dNEV
    dXYZ = x(1:3)' - TruePos;
    dNEV(i,:) = [gs, xyz2topo(dXYZ, TruePos_gd(1), TruePos_gd(2))];

    %% 결과 command 창에 표출
%     fprintf('Epoch = %d : delta x = %0.2f, delta y = %0.2f, delta z = %0.2f \n',...
%         i, TruePos(1) - x(1), TruePos(2) - x(2), TruePos(3) - x(3))
    fprintf('Epoch = %d : H error = %1.2fm: V error = %1.2fm : 3D error = %1.2fm \n',...
        i, norm(dNEV(i,2:3)), dNEV(i,4), norm(dNEV(i,2:4)))
    
    %% 산출값 저장
    % N 저장
    N_whole(i,:) = [gs, N'];
    % P 저장
    for p = 1: length(Pp(:,1))
        P_whole(i,1) = gs;
        P_whole(i,p+1) = P(p,p);
    end
    %% 실시간 PLOT
%     rt_Plot(FinalTTs, gs, x(1:3)', TruePos)

end


%% P plot
figure(22)
P_tHour = mod(P_whole(:,1), 86400); P_tHour = P_tHour/3600;
suptitle('P(x), P(y), P(z)')
subplot(3,1,1)
plot(P_tHour, P_whole(:,2),'b.:'); hold on; grid on;
axis([min(P_tHour) max(P_tHour) 0, 1]);
set(gca,'YTick',[0:0.2:1]);
ylabel('P_x');
sup2 = subplot(3,1,2);
plot(P_tHour, P_whole(:,3), 'b.:'); hold on; grid on;
axis([min(P_tHour) max(P_tHour) 0, 1]);
set(gca,'YTick',[0:0.2:1]);
ylabel('P_y');
sup3 = subplot(3,1,3);
plot(P_tHour, P_whole(:,4), 'b.:'); hold on; grid on; 
axis([min(P_tHour) max(P_tHour) 0, 1]);
set(gca,'YTick',[0:0.2:1]);
ylabel('P_z');
xlabel('Hours')

figure(33)
subplot(NoSats-1,1,NoSats-1)
xlabel('Hours')
for i = 1:NoSats-1
    prn = SatsList_OS(i);
    subplot(NoSats-1,1,i)
    plot(P_tHour, P_whole(:,4+i),'b.:'); hold on; grid on;
    xlim([min(P_tHour) max(P_tHour)])
    legend(['PRN : #',num2str(prn)])
end
        
%% N/C1 plot
for i=1:NoSats-1
    prn = SatsList_OS(i);
    QMRv_C1m = QMRv_C1(1:find(QMRv_C1(:,1) == FinalTTs(end)),:);
    C1_data = QMRv_C1(find(QMRv_C1m(:,2) == prn),:);
    N1_tHour = mod(N_whole(:,1), 86400); N1_tHour = N1_tHour/3600;
    C1_tHour = mod(C1_data(:,1), 86400); C1_tHour = C1_tHour/3600;
    figure(prn)
    title(['N vs C1(prn : ',num2str(prn),')']);
    yyaxis left
    plot(N1_tHour, N_whole(:,i+1),'b.'); hold on; grid on;
    xlim([min(C1_tHour) max(C1_tHour)]);
    xlabel('Hours')
    yyaxis right
    plot(C1_tHour, C1_data(:,4),'r.'); hold on; grid on;
    xlim([min(C1_tHour) max(C1_tHour)]);
    ylabel('Pseudorange(m)')
    legend('N','C1')
end

%% 오차 Plot
% [dXYZ, dNEV] = PosTErrorsJOON(estm(:,1), TruePos, estm(:,2:4));
[dXYZ, dNEV] = PosTErrorsRTK(estm(:,1), TruePos, estm(:,2:4), 0.3);