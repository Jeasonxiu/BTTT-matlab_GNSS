clear all;
close all;

%% �����Ŵ� �⺻ �Ķ���� ����
CCC = 299792458.;           	% ���� �ӵ� ��� ����. [m/s]
Freq_L1 = 1575.42e6;    Lambda = CCC/Freq_L1;    % GPS L1��ȣ ���ļ�/����
elecut = 15;

%% Base, Rover ����������
% ���� YY, DOY
% YY = '18'; DOY = '164';     % JAVAD-JAVAD
% YY = '18'; DOY = '159';     % JAVAD-JAVAD
YY = '18'; DOY = '165';     % JAVAD-JAVAD
% ���������� �ε�
% QMBs = load('QDSAJ_18164');    % 18164 JAVAD-A
% QMRv = load('QDSBS9_18164');    % 18164 S8
% QMBs = load('QDSAJ_18159');    % 18164 JAVAD-A
% QMRv = load('QDSBS1_18159');    % 18164 S8
% QMBs = load('QDAJA_18165_a');    % 18165 JAVAD-A
% QMRv = load('QDBS9_18165_a');    % 18165 S9(GEO++) L1 : cycle

% load('18165_rt_data.mat');
% QMBs = QDAN8_18165;    % 18165 JAVAD-A
% QMRv = QDBN8_18165;    % 18165 S9(GEO++) L1 : cycle
load('18169_rt_data.mat');
QMBs = QDAN8_18169;    % 18165 JAVAD-A
QMRv = QDBUB_18169;    % 18165 S9(GEO++) L1 : cycle
% QMRv = load('QDBUB_18169');    % 18165 S9(GEO++) L1 : cycle
% ����Ʈ�� gs �Ҽ��� ����
QMRv(:,1) = round(QMRv(:,1));
QMBs(:,1) = round(QMBs(:,1));
% 100���� GPS prn -> 100����
% QMBs(:,2) = QMBs(:,2)-100;
% QMRv(:,2) = QMRv(:,2)-100;    
QMBs_C1 = QMBs(QMBs(:,3) == 120, :);        % Base C1 ������ ����
QMBs_L1 = QMBs(QMBs(:,3) == 111, :);        % Base L1 ������ ����
QMRv_C1 = QMRv(QMRv(:,3) == 120, :);        % Base C1 ������ ����
QMRv_L1 = QMRv(QMRv(:,3) == 111, :);        % Base L1 ������ ����
QMRv_Pr = QMRv(QMRv(:,3) == 121, :);        % Base PrSigma ������ ����

%% L1 �Ÿ��� Cycle�� ��ȯ
% QMRv_L1(:,4) = QMRv_L1(:,4)/Lambda;       % GEO++������ Rinex������ ��Ȱ��ȭ

FinalTTs = intersect(unique(QMBs_L1(:,1)), unique(QMRv_L1(:,1)));
% FinalTT = FinalTTs(1:3300); FinalTTs = FinalTT;
%% �׹��޽��� �ε�
FileNav = strcat('brdc',DOY,'0.',YY,'n');   % GPS ��۱˵��� ���� �о���̱� =>  'brdc0000.00n'
eph = ReadEPH(FileNav); % �ش� brdc �����͸� �о����.
% load('eph180614.mat');
load('eph180618.mat');
eph(:,7) = sqrt(eph(:,7));
%% Base, Rover TruePos ����
% TruePos_Bs = [-3041235.578 4053941.677 3859881.013];        % A Point
TruePos_Bs = [-3041233.29541206 4053940.04178001 3859882.08156593];        % A Point
TruePos = [-3041241.741 4053944.143 3859873.640];           % B point
TruePos_gd = xyz2gd(TruePos);                               % TruePos geodetic coordinates
%% Kalman-filter �ʱⰪ ����
P = eye(3) * 500;                    % P xyz �ʱⰪ ����
%% Rover �ʱ���ǥ
AppPos = TruePos + 1;               % ������ǥ�� 1m ���� ����
AppPos = [-3041247.59543018,4053952.22498431,3859881.01792572];
% AppPos = [-3041247.16016828,4053947.71440359,3859872.00161267]; % 180515 S8 ��ü ù��ǥ
% �������� ���� �� �ι� ��ǥ ������ ���� x �� ����
x = AppPos';

%% GPS-RTK �������� ����
% �������� �ʱⰪ ����
Obs_Bs = zeros(32,3);
Obs_Rv = zeros(32,3);
SatsInfo_before = zeros(32,1);
SatsInfo_Bs = zeros(32,7);
SatsInfo_Rv = zeros(32,7);
N_before = zeros(32,1);
P_before = zeros(32,1);
x_before = zeros(32,1);
% �������� �ʱⰪ ����
idx_RS = 0;
prn_RS = 0;
for i = 1:length(FinalTTs)
%     for i = 1:2400
%         for i = 2217:2240
    vec_site = x(1:3)';
    gs = FinalTTs(i);                               % ���� Epoch gs
    % Base, Rover ���� Epoch ���������� ����
    QMBs_C1_1e = QMBs_C1(QMBs_C1(:,1) == gs, :);    % ���� gs�� C1 ������
    QMBs_L1_1e = QMBs_L1(QMBs_L1(:,1) == gs, :);    % ���� gs�� L1 ������
    QMRv_C1_1e = QMRv_C1(QMRv_C1(:,1) == gs, :);    % ���� gs�� C1 ������
    QMRv_L1_1e = QMRv_L1(QMRv_L1(:,1) == gs, :);    % ���� gs�� L1 ������
    QMRv_Pr_1e = QMRv_Pr(QMRv_Pr(:,1) == gs, :);    % ���� gs�� PrSigma ������
    Sats = intersect(QMBs_L1_1e(:,2),QMRv_L1_1e(:,2));
    SatsInfo_Bs = zeros(32,7);
    SatsInfo_Rv = zeros(32,7);
    % ���������� ���� ����ġ, ���� ��ǥ ��, ������, ���� ����
    if prn_RS ~= 0 & isempty(find(Sats(:) == prn_RS))
        idx_RS = 0;
        SatsInfo_before = zeros(32,1);
        N_before = zeros(32,1);
        P_before = zeros(32,1);
        x_before = zeros(32,1);
%         continue;
    end
    for j = 1:length(Sats)
        prn = Sats(j);

        % ���� gs�� prn�� �´� eph column ����
        icol = PickEPH(eph,prn,gs);
        toe = eph(icol, 8); a = eph(icol, 19); b = eph(icol, 20); c = eph(icol, 21); Tgd = eph(icol, 23);
        %% Base �������� �� ���� ���� ��� ����
        % Base �������� ��� ����
        C1_Bs = QMBs_C1_1e(QMBs_C1_1e(:,2) == prn,4);      % C1 observation value(Base)
        L1_Bs = QMBs_L1_1e(QMBs_L1_1e(:,2) == prn,4);      % L1 observation value(Base)
        Obs_Bs(prn,:) = [prn, C1_Bs, L1_Bs];
        % Base ���� ��ġ ���
        STT_Bs = GetSTTbrdc(gs, prn, eph, TruePos_Bs);
        tc_Bs = gs - STT_Bs;
        vec_sat_Bs = GetSatPosNC(eph, icol, tc_Bs);
        vec_sat_Bs = RotSatPos(vec_sat_Bs, STT_Bs);      % ���� ��ǥ
        vec_rho_Bs = vec_sat_Bs - TruePos_Bs';
        rho_Bs = norm(vec_rho_Bs);                    % Rho ����
        % Base ��ġ������ ���� ������, ���� ���
        Bs_gd = xyz2gd(TruePos_Bs);             % ������ ���� ����� ���� Base ��ǥ geodetic ��ȯ
        [az,el] = xyz2azel(vec_rho_Bs,Bs_gd(1),Bs_gd(2));
        if el > elecut
            SatsInfo_Bs(prn,:) = [prn, vec_rho_Bs', rho_Bs, az, el];
        end
        %% Rover �������� �� ���� ���� ��� ����
        % Rover �������� ��� ����
        C1_Rv = QMRv_C1_1e(QMRv_C1_1e(:,2) == prn,4);      % C1 observation value(Rover)
        L1_Rv = QMRv_L1_1e(QMRv_L1_1e(:,2) == prn,4);      % L1 observation value(Rover)
        Obs_Rv(prn,:) = [prn, C1_Rv, L1_Rv];
        % Rover ���� ��ġ ���
        STT = GetSTTbrdc(gs, prn, eph, vec_site);
        tc = gs - STT;
        vec_sat = GetSatPosNC(eph, icol, tc);
        vec_sat = RotSatPos(vec_sat, STT);      % ���� ��ǥ
        vec_rho = vec_sat - vec_site';
        rho = norm(vec_rho);                    % Rho ����
        % Rover ������ġ������ ���� ������, ���� ���
        Rv_gd = xyz2gd(vec_site);             % ������ ���� ����� ���� Rover ��ǥ geodetic ��ȯ
        [az,el] = xyz2azel(vec_rho,Rv_gd(1),Rv_gd(2));
        if el > elecut
            SatsInfo_Rv(prn,:) = [prn, vec_rho', rho, az, el];
        end
    end
    
    %% �������� ����
    if idx_RS == 0
        idx_RS = SatsInfo_Bs(find(SatsInfo_Bs(:,7) == max(SatsInfo_Bs(:,7))));
        prn_RS = SatsInfo_Bs(idx_RS,1); 
    end
    % �������� ����Ʈ�� ������ȣ�� �°� �籸��
    SatsList = intersect(SatsInfo_Bs(:,1), SatsInfo_Rv(:,1));
    SatsList = SatsList(find(SatsList > 0));
    k = 1;
    SatsInfo = [];
    for l=1:32
        prn = l;
        if k <= length(SatsList)
            if prn == SatsList(k) 
                SatsInfo(l,1) = prn;
                k = k + 1;
            end
        else
            SatsInfo(l,1) = 0;
        end
    end
    % ���������� ������ ������ ���� ����
    SatsList_OS = SatsList(SatsList(:) ~= prn_RS & SatsList(:) ~= 0);         % ���������� ������ �������
    %% ù epoch Ȥ�� ������ rise,set ��Ȳ(rise, set ������)
    if length(find((SatsInfo - SatsInfo_before) ~= 0))
        disp('rise/set');
        riseset(i,1) = gs;
    end
    % P, Q, A �ʱⰪ ���� 
    if length(find((SatsInfo - SatsInfo_before) ~= 0))
        P = P_initial(SatsInfo, SatsInfo_before, prn_RS, P(1:3,1:3), P_before);             % ù Epoch �� rise/set ���� P �ʱⰪ ����
    end
    Q = MakeQkt(length(SatsList_OS));
    A = eye(length(SatsList_OS)+3);                            % xyz + ���������� ������ NoSats
    % ��ȣ���� �ʱⰪ ����
    if length(find((SatsInfo - SatsInfo_before) ~= 0))
        N = N_initial(SatsInfo, SatsInfo_before, Obs_Bs, Obs_Rv, prn_RS, N_before);             % ù Epoch ���� N �ʱⰪ ����
    else
        N = N_before(find(N_before ~= 0));
        N_new = x(4:end);
    end
    % x ��迭
    if length(find((SatsInfo - SatsInfo_before) ~= 0))
        x = x_initial(SatsInfo, SatsInfo_before, prn_RS, x, x_before);             % ù Epoch ���� N �ʱⰪ ����
    end
    %% R: measurement noise
    R = MakeRrtk(SatsList_OS, SatsInfo_Rv);     % ���� ������ �̿��� R ��� ���� : ���ؼ�
%     R = makeR_KF(SatsList_OS, SatsInfo_Bs, SatsInfo_Rv, idx_RS, elecut);
%     R = MakeRrtk_Pr(SatsList_OS, SatsInfo_Rv,QMRv_Pr_1e);      % ������ PrSigmaM�� �̿��� R ��� ���� : ���ؼ�
%     R = MakeRkt(SatsInfo_Rv);     % ���� ������ �̿��� R ��� ���� : ������
% H matrix �ʱ�ȭ
    H=[]; y=[];
    %% �������� ��갪 ����
    for k=1:length(SatsList_OS)
        prn_OS = SatsList_OS(k);
        idx_OS = find(SatsInfo_Rv(:,1) == prn_OS);
        x(3+k) = N(k);
        if prn_OS == prn_RS
            continue;
        end
        % com ����
        com = (SatsInfo_Bs(idx_RS,5)-SatsInfo_Bs(idx_OS,5)) -...
            (SatsInfo_Rv(idx_RS,5)-SatsInfo_Rv(idx_OS,5)) +...
            Lambda*(N(k));
        COM(k,1) =com;
        obs = ((Obs_Bs(idx_RS,3)-Obs_Rv(idx_RS,3)) - (Obs_Bs(idx_OS,3)-Obs_Rv(idx_OS,3)))*Lambda;
        y(k,1) = obs -com;
        % H ��� ����
        H(k,1) = (SatsInfo_Rv(idx_RS,2)/SatsInfo_Rv(idx_RS,5)) - (SatsInfo_Rv(idx_OS,2)/SatsInfo_Rv(idx_OS,5));
        H(k,2) = (SatsInfo_Rv(idx_RS,3)/SatsInfo_Rv(idx_RS,5)) - (SatsInfo_Rv(idx_OS,3)/SatsInfo_Rv(idx_OS,5));
        H(k,3) = (SatsInfo_Rv(idx_RS,4)/SatsInfo_Rv(idx_RS,5)) - (SatsInfo_Rv(idx_OS,4)/SatsInfo_Rv(idx_OS,5));
        % ��ȣ���� ��̺а�
        H(k, 3+k) = Lambda;
    end

    %% Kalman-filter ����
    xp = A*x;
    Pp = A*P*A' + Q;
    K = Pp*H'/(H*Pp*H' + R);
    % ��ǥ ����
    x = xp + K*(y);
    % ���л� ��� ����
    P = Pp - K*H*Pp;
    estm(i,:) = [gs, x(1:3)'];
    
    %% N, P, x ��迭 �� ���� epoch �������� ����
    N_before = makeNbefore(x(4:end), prn_RS, SatsInfo);
    P_before = makePbefore(P, prn_RS, SatsInfo);
    x_before = makexbefore(x, prn_RS, SatsInfo);
    SatsInfo_before = SatsInfo;
    %% dXYZ, dNEV
    dXYZ = x(1:3)' - TruePos;
    dNEV(i,:) = [gs, xyz2topo(dXYZ, TruePos_gd(1), TruePos_gd(2))];

    %% ��� command â�� ǥ��
    NoSats = length(SatsInfo(find(SatsInfo ~= 0)));
    
%     fprintf('Epoch = %d : delta x = %0.2f, delta y = %0.2f, delta z = %0.2f \n',...
%         i, TruePos(1) - x(1), TruePos(2) - x(2), TruePos(3) - x(3))
    fprintf('Used Sats = %d : Epoch = %d : H error = %1.2fm: V error = %1.2fm : 3D error = %1.2fm \n',...
        NoSats, i, norm(dNEV(i,2:3)), dNEV(i,4), norm(dNEV(i,2:4)))
%     fprintf('Used Sats = %d : gs = %d : H error = %1.2fm: V error = %1.2fm : 3D error = %1.2fm \n',...
%         NoSats, gs, norm(dNEV(i,2:3)), dNEV(i,4), norm(dNEV(i,2:4)))
    List = SatsInfo(find(SatsInfo ~= 0))';
%     disp(['Used Sats = ', num2str(NoSats)]);
%     disp(['Used Sats = ', num2str(List)]);
    %% ���Ⱚ ����
    % N ����
%     N_whole(i,:) = [gs, N'];
    % P ����
    for p = 1: length(Pp(:,1))
        P_whole(i,1) = gs;
        P_whole(i,p+1) = P(p,p);
    end
    %% �ǽð� PLOT
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

% figure(33)
% subplot(NoSats-1,1,NoSats-1)
% xlabel('Hours')
% for i = 1:NoSats-1
%     prn = SatsList_OS(i);
%     subplot(NoSats-1,1,i)
%     plot(P_tHour, P_whole(:,4+i),'b.:'); hold on; grid on;
%     xlim([min(P_tHour) max(P_tHour)])
%     legend(['PRN : #',num2str(prn)])
% end
%         
% %% N/C1 plot
% for i=1:NoSats-1
%     prn = SatsList_OS(i);
%     QMRv_C1m = QMRv_C1(1:find(QMRv_C1(:,1) == FinalTTs(end)),:);
%     C1_data = QMRv_C1(find(QMRv_C1m(:,2) == prn),:);
%     N1_tHour = mod(N_whole(:,1), 86400); N1_tHour = N1_tHour/3600;
%     C1_tHour = mod(C1_data(:,1), 86400); C1_tHour = C1_tHour/3600;
%     figure(prn)
%     title(['N vs C1(prn : ',num2str(prn),')']);
%     yyaxis left
%     plot(N1_tHour, N_whole(:,i+1),'b.'); hold on; grid on;
%     xlim([min(C1_tHour) max(C1_tHour)]);
%     xlabel('Hours')
%     yyaxis right
%     plot(C1_tHour, C1_data(:,4),'r.'); hold on; grid on;
%     xlim([min(C1_tHour) max(C1_tHour)]);
%     ylabel('Pseudorange(m)')
%     legend('N','C1')
% end

%% ���� Plot
% [dXYZ, dNEV] = PosTErrorsJOON(estm(:,1), TruePos, estm(:,2:4));
[dXYZ, dNEV] = PosTErrorsRTK(estm(:,1), TruePos, estm(:,2:4), 0.3);