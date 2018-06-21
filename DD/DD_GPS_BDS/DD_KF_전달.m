%% Double Difference Algorithm using Carrier Phase (Single Frequency)
clear all;close all;%clc;

%% ����� ������ ����
Maximum_Iteration = 1;    % �ִ� �ݺ�Ƚ�� ������
stopConditon = 1e-5;      % Iteration �Ӱ谪
eleCut = 15;                % �Ӱ���� ����

%% ������ ����� QM file ����.
FileQM_Bs ='QGSDA_17025H1.qm';   %  GSDA (���굿 ���� A)
FileQM_Rv ='QGSDB_17025H1.qm';   %  GSDB (���굿 ���� B)
%% Part 1. Between-Satellite Single Difference Positioning algorithms.

CCC = 299792458.;           % ���� �ӵ� ��� ����. [m/s]
Freq_L1 = 1575.42e6;        % L1 ���ļ�
Lambda_L1 = CCC/Freq_L1;    % L1 ����

[arrQM_Bs, FinalPRNs_Bs, FinalTTs_Bs] = ReadQM(FileQM_Bs);  % ������ ����� QM file�� �о����.
[arrQM_Rv, FinalPRNs_Rv, FinalTTs_Rv] = ReadQM(FileQM_Rv);  % ������ ����� QM file�� �о����.
FinalTTs = unique(intersect(FinalTTs_Bs,FinalTTs_Rv));      % �� ���� �������� ���� ���Žð� ����

ObsType = 111;                                  % ����ġ Ÿ�� :  L1 (CP)
ObsType_CA = 120; % �ð����� �������� �ڵ� ����ġ ��� - 120: CA

ObsType_Selected_QM_Bs = SelectQM(arrQM_Bs, ObsType);       % ������ ����� ����ġ�� ����. (in QM_Bs)
ObsType_Selected_QM_Rv = SelectQM(arrQM_Rv, ObsType);       % ������ ����� ����ġ�� ����. (in QM_Rv)

% �ش� QM file�� �����ϴ� brdc �����͸� ����.
YY = FileQM_Bs(7:8);
DOY = FileQM_Bs(9:11);
FileNav = strcat('brdc',DOY,'0.',YY,'n');

eph = ReadEPH(FileNav); % �ش� brdc �����͸� �о����.

[gw, gd] = ydoy2gwgd(str2num(YY), str2num(DOY)); % ��꿡 �̿��� GPS week number ����.

% ����� ��ǥ
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

%% Kalman Filter�� �̿��� ���� ����
% ������ ���� ������ ����

NoEpochs = length(FinalTTs); % Number of Epochs.

Estimated_Position = zeros(NoEpochs, 5);   % ����ġ ������ ���� ���. (�� ������� �ð�, X, Y, Z, �ð����)
No_Est = 0;

tic
%% Į������
    Sats_for_Process = FinalPRNs_Bs;
    NoSats = length(Sats_for_Process);
    
    %%
    x = zeros(3+(NoSats-1),1); % ���� ��� ����� ���� �ʱⰪ ����.
    x(1:3) = AppPos_Rv';
        
    A = eye(3+(NoSats-1));
    A = 1 .* A;

    Q = eye(3+(NoSats-1));
    P = eye(3+(NoSats-1));

    save_P = zeros(1,(length(P))+2);
%% ����
for i = 1:NoEpochs
    %%
    IndexS_QM_Bs = find(ObsType_Selected_QM_Bs(:,1) == FinalTTs(i));   % �� �ð��� L1 measurements ���� ������ ���� index.
    IndexS_QM_Rv = find(ObsType_Selected_QM_Rv(:,1) == FinalTTs(i));   % �� �ð��� L1 measurements ���� ������ ���� index.
    
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
    
    gs = FinalTTs(i);               % ��꿡 �̿��� GPS second ����.
    
    rcvPos = AppPos_Rv;                % ���ű� xyz ��ǥ
    
    %% ��������(RS) ���� //  SatsAZEI - c1(gs), c2(prn), c3(az), c4(el)
    [SatsAZEl, indxRS]=pickRSel(gs,Sats_for_Process, eph, rcvPos, Correction);      % RS : Reference Satellite ��������
    RS_PRN = Sats_for_Process(indxRS);
    
    %% ��������(RS)�� ����ġ ����
    obs_RS_Bs = sub_QM_Bs(find(sub_QM_Bs(:,2)==RS_PRN),4);
    obs_RS_Rv = sub_QM_Rv(find(sub_QM_Rv(:,2)==RS_PRN),4);
    
    %% �ݺ� ��� ����
        
        % �ּ��ڽ¹� ����� ���� Matrix ����.
        H = zeros(1,3+(NoSats-1));
        y = zeros(1,1);             % y(����ġ-��갪) ����
        
        R = eye(NoSats-1,NoSats-1);  % R matrix ����
        
        % �ݺ���궧���� ���ο� x���� ���
        vec_site = x(1:3)';              % ���ű� xyz ��ǥ (�ݺ������ ������ ������ ���ŵ�.)
        
        %% �������� ��ǥ ���� ���
        % �ش� �ð��� ������ �����͸� �̿��ϱ� ���� column ��ġ ����.
        indxEPH_RS = PickEPH(eph, RS_PRN, gs);
        if indxEPH_RS == 0
            continue;
        end
        
        STT = obs_RS_Bs(i)/CCC;

        tc = gs - STT;
        vec_Sat_RS = GetSatPos_G(eph, indxEPH_RS, tc);      % �ش� ������ ��ġ ����
        
        vec_Sat_RS = RotSatPos(vec_Sat_RS, STT);
        
        vec_rho_RS_Bs = vec_Sat_RS - TruePos_Bs; 
        vec_rho_RS_Rv = vec_Sat_RS - vec_site; 
        rho_RS_Bs = norm(vec_rho_RS_Bs);
        rho_RS_Rv = norm(vec_rho_RS_Rv);

        No_OS = 0;
        
        %% �� ������ ���� ����ġ, ���ġ, H��İ��
        for j = 1:NoSats

            el_OS = SatsAZEl(j,4);

            OS_PRN = Sats_for_Process(j);         	  % �ش� ������ PRN
            
            if OS_PRN == RS_PRN
                continue
            else
                No_OS = No_OS+1;
            end
            
            %% DD����ġ ������Ʈ
            obs_OS_Bs = sub_QM_Bs(find((sub_QM_Bs(:,1)==gs)&(sub_QM_Bs(:,2)==OS_PRN)),4);    % y ��� ����� ���� �ش� �ð�, �ش� ������ ������(obs) ����.
            obs_OS_Rv = sub_QM_Rv(find((sub_QM_Rv(:,1)==gs)&(sub_QM_Rv(:,2)==OS_PRN)),4);    % y ��� ����� ���� �ش� �ð�, �ش� ������ ������(obs) ����.

            % �ش� �ð��� ������ �����͸� �̿��ϱ� ���� column ��ġ ����.
            indxEPH_OS = PickEPH(eph, OS_PRN, gs);   
            if indxEPH_OS == 0
                continue;
            end

            STT = obs_OS_Bs/CCC;

            tc = gs - STT;
            vec_Sat_OS = GetSatPos_G(eph, indxEPH_OS, tc);       % �ش� ������ ��ġ ����
            
            vec_Sat_OS = RotSatPos(vec_Sat_OS, STT);
        
            %% SD ���ġ ������Ʈ-�Ÿ� ���ġ�� ���� ����� ���� DD���ġ ���

            vec_rho_OS_Bs = vec_Sat_OS - TruePos_Bs;     % ������ ���űⰣ�� rho vector ����.
            vec_rho_OS_Rv = vec_Sat_OS - vec_site;     % ������ ���űⰣ�� rho vector ����.
            rho_OS_Bs = norm(vec_rho_OS_Bs);              % ��� ����� ���� �������κ��� ���ű������ �Ÿ� ���.
            rho_OS_Rv = norm(vec_rho_OS_Rv);              % ��� ����� ���� �������κ��� ���ű������ �Ÿ� ���.

            com = ((rho_RS_Bs - rho_OS_Bs) - (rho_RS_Rv - rho_OS_Rv));
            obs = ((obs_RS_Bs - obs_OS_Bs) - (obs_RS_Rv - obs_OS_Rv));
            
            if el_OS > eleCut %% - ������ �ٲ㰡�鼭 �غ���
                R(No_OS,No_OS) = 1;
            else
                R(No_OS,No_OS) = ( 1/sind(el_OS) );
            end
                
            y(No_OS,1) = obs - com;    % ����ġ - ��갪
            
            % H Matrix ����.
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
        Estimated_Position(No_Est,1) = gs;            % ����ġ ����. (�� ������� �ð�, X, Y, Z, �ð����)
        Estimated_Position(No_Est,2:5) = x(1:4);      % ����ġ ����. (�� ������� �ð�, X, Y, Z, �ð����)
 
        for ppppp = 1:(length(P))
            save_P(No_Est,ppppp) = P(ppppp,ppppp);
        end
        save_P(No_Est,((length(P))+1)) = sub_QM_Bs(1,1);
        save_P(No_Est,((length(P))+2)) = NoSats;

end
toc


%% Part 2. Positioning Accuracy Plotting 
Estimated_Position = Estimated_Position(1:No_Est, :);
estmTT = Estimated_Position(:, 1);        % ���� ����� �����ϴ� Time Tags (gps second)
estmPos = Estimated_Position(:, 2:4);     % ���� ���� ���ű���ġ
NoPos = length(estmPos);    % ���� ��� �������� ����.

geodetic_true = xyz2gd(TruePos_Rv); % TRUE POSITION Latitude, Longitude
TrueLat = geodetic_true(1);      % TRUE POSITION Latitude
TrueLon = geodetic_true(2);      % TRUE POSITION Longitude

% delta X, delta Y, delta Z ����
dXYZ = zeros(NoPos,3);
for k = 1:NoPos
    dXYZ(k,:) = estmPos(k,:) - TruePos_Rv;
end

dNEV = xyz2topo(dXYZ, TrueLat, TrueLon); % delta XYZ�� delta NEV�� ��ȯ


% -------------------------------------------------------------------------
% RMSE ����
dN = dNEV(:,1); 
dE = dNEV(:,2); 
dNE = sqrt(dN.^2 + dE.^2);       
rmsH = sqrt(mean(dNE.^2));  % ���� ���� ����
Last_H = sqrt(dNE(end)^2);  % ���� ���� ����

mean_dN = sqrt(mean(dN.^2));            % ���Ϲ��� ���� ���� ����
Last_dN = sqrt(dN(end)^2);            % ���Ϲ��� ���� ���� ����
mean_dE = sqrt(mean(dE.^2));            % �������� ���� ���� ����
Last_dE = sqrt(dE(end)^2);            % �������� ���� ���� ����

dV = dNEV(:,3);
rmsV = sqrt(mean(dV.^2));   % ���� ���� ����
Last_V = sqrt(dV(end)^2);   % ���� ���� ����

d3 = sqrt(dN.^2 + dE.^2 + dV.^2); 
rms3 = sqrt(mean(d3.^2));   % 3���� ���� ����
Last_3 = sqrt(d3(end)^2);   % 3���� ���� ����


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
    tHour = tHour/3600;          % ��� �ð������� Hour�� ����.
    
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
       
    
