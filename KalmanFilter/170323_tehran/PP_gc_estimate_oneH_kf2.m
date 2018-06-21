
clear all; close all;

% % %
% % QMfile = 'QM170125_A';
% % QMfile = 'T4_170323_k500';
% % QMfile = 'T4_GC';
% QMfile = 'QM170323_Bs';        % ���� bs
% % QMfile = 'QMfile';
% % navfile = 'brdm0250.17p';
% % navfile = 'hour0450.17n';
% navfile = 'brdm0820.17p';
% % TruePos = [-3041235.578 4053941.677 3859881.013];
% TruePos = [-3052725.48122377,4042747.35562486,3862428.65171442];
% YY = 17; DOY =082;
% % TruePos = [-3058799.61420451,4083265.35912516,3814946.87192938];
%
% %% �Һ� ���� ����: ���� �ӵ�, ����ġ
% CCC = 299792458.;   % CCC = Speed of Light [m/s]
% %% �Ӱ���� ����
% eleCut = 15;
% %% ���� �ڵ鸵
% LeapSecBDS = 14;
%
% %% ����� ����ġ ����
% g_ObsType = 120; % gps C1
% g_ObsType_snr = 141;
% g_ObsType_dop = 131;
% c_ObsType = 220; % bds C1
% c_ObsType_snr = 241;
% c_ObsType_dop = 231;
%
% %% QM ���� �о�鿩�� ��ķ� �����ϰ�, ����� ����ġ ����
% [arrQM, FinalPRNs, FinalTTs] = ReadQM(QMfile);
% arrQM = arrQM(find(arrQM(:,3) < 300),:);
% QM = SelectQM_gc(arrQM, g_ObsType, c_ObsType);
% QM_snr = SelectQM_gc(arrQM, g_ObsType_snr, c_ObsType_snr);
% QM_dop = SelectQM_gc(arrQM, g_ObsType_dop, c_ObsType_dop);
% FinalTTs = unique(QM(:,1));
% % FinalTTs = FinalTTs(find(FinalTTs(:,1) > 1000),:);
%
% % %% eph ����
% [eph, trashPRN, trashT]=ReadEPH_all(navfile);
%
% %% �׹��޽����� �о�鿩�� ��ķ� �����ϰ�, Klobuchar �� ����
% gps_nav = strcat('brdm', num2str(DOY,'%03d'), '0.', num2str(YY,'%02d'), 'p');   %: Navigation RINEX file
% fid = fopen(gps_nav,'r');
% if fid == -1
%     al = zeros(4,1); be = zeros(4,1);
% else
%     [al, be] = GetALBE(gps_nav);
% end
%
% %% �׹��޽��� ���� �̸��� �̿��� YY, DOY ����
% % YY = str2num(navfile(length(navfile)-2:length(navfile)-1));
% % DOY = str2num(navfile(length(navfile)-7:length(navfile)-5));
% [gw, gd] = ydoy2gwgd(YY, DOY); %: GPS WEEK ����
%
% %% �ʱ���ǥ ȹ��
% AppPos = TruePos;
% gd = xyz2gd(AppPos); AppLat = gd(1); AppLon = gd(2);
%
% %% ������ �ʿ��� �ʱ�ġ ����
% Maxiter = 10;
% EpsStop = 1e-5;
% ctr = 1; ctr2 = 1;
% x = [AppPos ctr ctr2]; x = x';
%
% %% �������� ����
% NoEpochs = length(FinalTTs);
% estm = zeros(NoEpochs, 6);
% Scount = zeros(NoEpochs, 1);
load('PP_170323_gc.mat');
nEst = 0;

% % P �ʱ�ȭ
P = eye(5);
P(1:3,1:3) = P(1:3,1:3)*1;
P(4,4) = P(4,4)*0.5; P(5,5) = P(5,5)*0.5;

% Q �ʱ�ȭ
Q = eye(5);
Q(1,1) = 5;
Q(2,2) = 5;
Q(3,3) = 5;
Q(4,4) = 10000; Q(5,5) = 10000;

order = 2;

%useful values to index matrices
o1 = order;
o2 = order*2;
o3 = order*3;
sigmaq0 = 1;
sigmaq_vE = 100^2;
sigmaq_vN = 100^2;
sigmaq_vU = 200^2;
T0 = eye(o1) + diag(ones(o1-1,1),1)*1;
Z_1_om = zeros(1,o1-1);
I = eye(o3);
Z_o1_o1 = zeros(o1);
T = [T0      Z_o1_o1 Z_o1_o1;
    Z_o1_o1 T0      Z_o1_o1;
    Z_o1_o1 Z_o1_o1 T0];
% xkf = [x(1);0;x(2);0;x(3);0];
% x= [-3108688.89600000,4078504.77600000,3779792.40400000]';
% x= [-3108697.20600000,4078501.36200000,3779789.09700000]';
xkf = [x(1);0;x(2);0;x(3);0];
eleCut =5;

for j = 1:NoEpochs
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
            STT = GetSTTbrdm(gs, eph, icol, x(1:3)); % ��ȣ���޽ð� ���
            if prn > 100 && prn < 200
                tc = gs - STT ;             % ��ȣ���޽ð� ����
            elseif prn > 200 && prn < 300
                tc = gs - STT- LeapSecBDS;  % ��ȣ���� �ð� ����... bds���� �ݿ�
            end
            SatPos = GetSatPosNC_GC(eph, icol, tc); % ������ġ ����
            SatPos = RotSatPos(SatPos, STT); % ��������ȿ�� ���
            vec_rho = SatPos - vec_site';
            rho = norm(vec_rho);
            [az,el] = xyz2azel(vec_rho, AppLat, AppLon);
            if prn > 100 && prn < 200
                %                     dIono = 0;
                %                     dTrop = 0;
                dIono = klobuchar(al, be, gs, az, el, x(1:3)); % �̿��� ����(Klobuchar ��)
                dTrop = ZHD2SHD(gw, gs, TruePos, el, ZHD); % ����� ����
            elseif prn > 200 && prn < 300
                %                     dIono = 0;
                %                     dTrop = 0;
                dIono = klobuchar(al, be, gs, az, el, x(1:3)); % �̿��� ����(Klobuchar ��)
                dTrop = ZHD2SHD(gw, gs, TruePos, el, ZHD); % ����� ����
            end
            %             dRel=0;
            dRel = GetRelBRDC(eph, icol, tc); % ��뼺ȿ��
            dtSat = a + b*(tc - toe) + c*(tc - toe)^2 - Tgd + dRel; % �����ð���� ���... �׷������, ��뼺ȿ�� ����
            if el >=eleCut % �Ӱ����
                W = 1;
                %                 W = MakeW_elsnr(el,snr);
                
                elsnr(cnt2,1:5) = [gs, prn, el, snr, dop];
                
                if prn > 100 && prn < 200
                    %                         PRC = PickPRC(RTCM,prn,gs);
                    PRC = 0;
                    com = rho + x(4) - CCC*dtSat + dTrop + dIono - PRC; % gps ��갪
                    Com = rho - CCC*dtSat + dTrop + dIono - PRC; % gps ��갪
                    H(cnt2,:) = [-vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 1 0];
                elseif prn > 200
                    %                         PRC = PickPRC(RTCM,prn,gs);
                    PRC = 0;
                    com = rho + x(5) - CCC*dtSat + dTrop + dIono - PRC; % bds ��갪
                    Com = rho - CCC*dtSat + dTrop + dIono - PRC; % gps ��갪
                    H(cnt2,:) = [-vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 1 1];
                end
                y = obs - com;
                Y0(cnt2,1) = obs;
                COM(cnt2,1) = Com;
                Y(cnt2,1) = y;
                cnt2 = cnt2 + 1;
                NoSatsUsed = NoSatsUsed + 1;
            end
        end
        velocity = GetVelDop(gs, elsnr(:,2), elsnr(:,5), eph, x(1:3)');
        W = PP_cofactor_matrix(elsnr(:,3), elsnr(:,4), 4);
        HTH = H(1:cnt2-1,:)'*inv(W(1:cnt2-1,1:cnt2-1))*H(1:cnt2-1,:);
        HTy = H(1:cnt2-1,:)'*inv(W(1:cnt2-1,1:cnt2-1))*Y(1:cnt2-1,:);
        xhat = inv(HTH) * HTy;
        y0 = x(1:3) + xhat(1:3);
        x(1:3)= y0;
        y_hat = H*xhat + COM;
        v_hat = Y0 - y_hat;
        v_hat = v_hat*0.1;
        invQ=diag((diag(W).^-1));
        sigma02_hat(1,1) = (v_hat'*(W^-1)*v_hat) / (cnt2-1-3);
        sigma02_hat(1,2) = (v_hat'*(invQ)*v_hat);
        sigma02_hat(1,3) = (cnt2-1-5);
%         sigma02_hat(1,1) = sigma02_hat(1,2)/sigma02_hat(1,3); 
        Cxx = sigma02_hat(1) * (HTH^-1);
        cov_XR  = Cxx(1:3,1:3);
        sigma2_XR = diag(cov_XR);
        if j == 1
            Cee(:,:) = zeros(o3);
            Cee(1,1) = sigma2_XR(1);
            Cee(o1+1,o1+1) = sigma2_XR(2);
            Cee(o2+1,o2+1) = sigma2_XR(3);
            Cee(2:o1,2:o1) = sigmaq0 * eye(o1-1);
            Cee(o1+2:o2,o1+2:o2) = sigmaq0 * eye(o1-1);
            Cee(o2+2:o3,o2+2:o3) = sigmaq0 * eye(o1-1);
        end
        Cvv(o1,o1) = sigmaq_vE;
        Cvv(o2,o2) = sigmaq_vN;
        Cvv(o3,o3) = sigmaq_vU;
        Cvv([o1 o2 o3],[o1 o2 o3]) = local2globalCov(Cvv([o1 o2 o3],[o1 o2 o3]), xkf([1 o1+1 o2+1]));
        %             Cvv([o1 o2 o3],[o1 o2 o3]) = local2globalCov(Cvv([o1 o2 o3],[o1 o2 o3]), x);
        Cnn =cov_XR(1:3,1:3);
        P = Cee;
        Q = Cvv;
        Hh = [1 Z_1_om 0 Z_1_om 0 Z_1_om;
            0 Z_1_om 1 Z_1_om 0 Z_1_om;
            0 Z_1_om 0 Z_1_om 1 Z_1_om];
        
        K = T*Cee*T' + Cvv;
        G = K*Hh' * (Hh*K*Hh' + Cnn)^-1;
        Xhat = (I-G*Hh)*xkf + G*y0;
        xkf = T * Xhat;
        Cee = (I-G*Hh)*K;
        nEst = nEst + 1;
        estm(nEst,1) = gs;
        estm(nEst,2:4) = xkf(1:2:5);
        estm(nEst,7) = NoSatsUsed;
        Scount(nEst,1) = NoSatsUsed;
        fprintf('%8d : %3d : %8.2f : %8.2f\n', j, iter, xkf(1)' - TruePos(1), xkf(3)' - TruePos(2));
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

figure()
hold on; grid on;
% plot(Base_gd(:,3), Base_gd(:,2), 'b.')
plot(estm_gd(1:2500,3), estm_gd(1:2500,2), 'r.','Markersize',10)
% plot(estm_gd(:,3), estm_gd(:,2), 'r.','Markersize',10)
plot_google_map