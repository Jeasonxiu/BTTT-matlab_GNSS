%% GPS & BDS ���� ���� �ڵ� (DGNSS)
% tic;
clear all; close all;
warning off;
%% ��� ����
CCC = 299792458.;   % CCC = Speed of Light [m/s]
%% QM ����, ��¥ �ڵ鸵

%-----------170116 �׽�Ʈ---------%
RTCM = load('PPS1_170302.t41');
DOY = 061;YY  = 17;
%% Navigation file load
FileNav = strcat('brdm', num2str(DOY,'%03d'), '0.', num2str(YY,'%02d'), 'p');   %: Navigation RINEX file
[eph, trashPRN, trashT]=ReadEPH_all(FileNav);
% load('PTH_DGNSS.mat');

for LR=15:16
    % for LR=1:24
    LR
    tic;
    switch LR
        %----���� Type A----%
        case 1
            FileQM='PTH2_170302_sedan_left1';
        case 2
            FileQM='PTH2_170302_sedan_right1';
            %----���� Type B----%
        case 3
            FileQM='PTH2_170302_sedan_left2';
        case 4
            FileQM='PTH2_170302_sedan_right2';
        case 5
            FileQM='PTH2_170302_sedan_left3';
        case 6
            FileQM='PTH2_170302_sedan_right3';
            %----SUV Type A----%
        case 7
            FileQM='PTH2_170302_SUV_left1';
        case 8
            FileQM='PTH2_170302_SUV_right1';
            %----SUV Type B----%
        case 9
            FileQM='PTH2_170302_SUV_left2';
        case 10
            FileQM='PTH2_170302_SUV_right2';
        case 11
            FileQM='PTH2_170302_SUV_left3';
        case 12
            FileQM='PTH2_170302_SUV_right3';
        case 13
            FileQM='PTH2_170302_LINE_left';
        case 14
            FileQM='PTH2_170302_LINE_right';
        case 15
            FileQM='PTH2_170302_static_left';
        case 16
            FileQM='PTH2_170302_static_right';
    end
    %-----------------------------------------------------%
    [arrQM, FinalPRNs, FinalTTs] = ReadQM(FileQM);
    
    %% �̻� ���� ����
%     switch LR
%         
%         case 9
%             idx_wrong = [min(find(arrQM(:,1) == 110332)),max(find(arrQM(:,1) == 110362))];      % test4_1_L
%             arrQM_1st = arrQM(1:idx_wrong(1)-1,:);
%             arrQM_2nd = arrQM(idx_wrong(2)+1:end,:);
%             arrQM_10 = arrQM(idx_wrong(1):idx_wrong(2),:);
%             arrQM_10 = arrQM_10(find(arrQM_10(:,2) ~= 25),:);       % test4_1_R
%             arrQM = [arrQM_1st; arrQM_10; arrQM_2nd];
%         case 10
%             idx_wrong = [min(find(arrQM(:,1) == 363228)),max(find(arrQM(:,1) == 363247))];      % test4_1_R
%             arrQM_1st = arrQM(1:idx_wrong(1)-1,:);
%             arrQM_2nd = arrQM(idx_wrong(2)+1:end,:);
%             arrQM_10 = arrQM(idx_wrong(1):idx_wrong(2),:);
%             arrQM_10 = arrQM_10(find(arrQM_10(:,2) ~= 25),:);       % test4_1_R
%             arrQM = [arrQM_1st; arrQM_10; arrQM_2nd];
%     end
    % break
    yyS = num2str(YY,'%02d');
    doyS = num2str(DOY,'%03d');
    [gw, gd] = ydoy2gwgd(YY, DOY); %: GPS WEEK ����
    [bw, bd] = ydoy2bwbd(YY, DOY); %: BDS WEEK ����
    
    %% ��ǥ ����
%     TruePos = [-3058799.61420451,4083265.35912516,3814946.87192938]; %170116 ������ �ʱ���ǥ
    TruePos = [-3058726.836  4083297.744  3814946.173];
    
    %% ���� �ڵ鸵
    LeapSecBDS = 14;
    
    %% brdm���Ͽ��� al/be�� ���� ������ brdc���Ͽ��� ������
    al=[0 0 0 0];   be=[0 0 0 0];
    
    %% ����� ����ġ ����
    g_ObsType = 120; % gps C1
    g_ObsType_snr = 141;
    c_ObsType = 220; % bds C1
    c_ObsType_snr = 241;
    
    %% QM �غ�
    QM = SelectQM_gc(arrQM, g_ObsType, c_ObsType);
    QM_snr = SelectQM_gc(arrQM, g_ObsType_snr, c_ObsType_snr);
    FinalPRNs = unique(QM(:,2));
    FinalTTs = unique(QM(:,1));
    FinalTTs = FinalTTs(find(FinalTTs(:,1) > 1000),:);
    % break
    
    %% ��ǥ �ʱ�ġ ���� �� ���浵 ��ȯ
    AppPos = TruePos;
    gd = xyz2gd(AppPos); AppLat = gd(1); AppLon = gd(2);
    
    %% ������ ���� �Ű����� ����
    Maxiter = 10;
    EpsStop = 1e-4;
    ctr = 1; ctr2 = 1;
    eleCut = 15;
    x = [AppPos ctr ctr2]';
    %%
    NoEpochs = length(FinalTTs);
    estm = zeros(NoEpochs, 6);
    Scount = zeros(NoEpochs, 1);
    nEst = 0;
    SNR = 0;
    cnt=1;
    for j = 1:NoEpochs % 60:65%
        x = [AppPos ctr ctr2]';
        for iter = 1:Maxiter
            gs = FinalTTs(j);
            HTH = zeros(5,5);
            HTy = zeros(5,1);
            indexQM = find(QM(:,1) == gs);
            QM1e = QM(indexQM,:);
            QM1e_snr = QM_snr(indexQM,:);
            NoSats = length(QM1e);
            NoSatsUsed = 0;
            
            if NoSats < 6
                break
            end
            vec_site = x(1:3)';
            ZHD = TropGPTh(TruePos, gw, gs); %: TROP: GPT
            
            for i = 1:NoSats
                prn = QM1e(i,2);
                obs = QM1e(i,4);
                snr = QM1e_snr(i,4);
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
                    dIono = 0;
                    dTrop = 0;
                    %             dIono = ionoKlob(al, be, gs, az, el, x(1:3)); % �̿��� ����(Klobuchar ��)
                    %             dTrop = ZHD2SHD(gw, gs, TruePos, el, ZHD); % ����� ����
                elseif prn > 200 && prn < 300
                    dIono = 0;
                    dTrop = 0;
                    %             dIono = ionoKlob(al, be, gs, az, el, x(1:3)); % �̿��� ����(Klobuchar ��)
                    %             dTrop = ZHD2SHD(gw, gs, TruePos, el, ZHD); % ����� ����
                end
                %             dRel=0;
                dRel = GetRelBRDC(eph, icol, tc); % ��뼺ȿ��
                dtSat = a + b*(tc - toe) + c*(tc - toe)^2 - Tgd + dRel; % �����ð���� ���... �׷������, ��뼺ȿ�� ����
                if el >=eleCut % �Ӱ����
                    AZEL(cnt,1:4) = [prn, el, az, snr];
                    cnt=cnt+1;
                    %                 W = 1;
                    W = MakeW_elsnr(el,snr);
                    if iter == 1
                        SNR = SNR + snr;
                    end
                    if prn > 100 && prn < 200
%                         PRC = PickPRC(RTCM,prn,gs);
                                            PRC = 0;
                        com = rho + x(4) - CCC*dtSat + dTrop + dIono - PRC; % gps ��갪
                        H = [-vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 1 0];
                    elseif prn > 200
%                         PRC = PickPRC(RTCM,prn,gs);
                                                                    PRC = 0;
                        com = rho + x(5) - CCC*dtSat + dTrop + dIono - PRC; % bds ��갪
                        H = [-vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 0 1];
                    end
                    y = obs - com;
                    HTH = HTH + H'*W*H;
                    HTy = HTy + H'*W*y;
                    NoSatsUsed = NoSatsUsed + 1;
                end
            end
            P = inv(HTH);
            xhat = inv(HTH) * HTy;
            x = x + xhat;
            if norm(xhat) < EpsStop;
                nEst = nEst + 1;
                estm(nEst,1) = gs;
                estm(nEst,2:6) = x(1:5);
                estm(nEst,7) = NoSatsUsed;
                estm(nEst,8) = SNR;
                Scount(nEst,1) = NoSatsUsed;
                SNR = 0;
%                 fprintf('%8d : %8.2f \n', j, norm(x(1:3)'-TruePos));
                break;
            end
        end
    end
    estm = estm(1:nEst,:);
    Scount=Scount(1:nEst, :);
    
    
    fixed_tmp = ([FileQM(13:end),'=estm;']);
    eval(fixed_tmp);
    toc;
end

% figure(101);
% [dXYZ, dNEV]=PosErrors(estm(:,1), TruePos,estm(:,2:5),Scount);
% break
%
% %% �Է¹��� �ڷ��� ũ�� �� ������ �ش��ϴ� ���浵 ����
% NoPos = length(estm(:,1));
% gd = xyz2gd(TruePos); TrueLat = gd(1); TrueLon = gd(2);
% %% ������ XYZ�� ���� XYZ�� ���̰� ����
% dXYZ = zeros(NoPos,3);
% for k = 1:NoPos
%     gd(k,:) = xyz2gd(estm(k,2:4));
% end
% estm_L=estm;
% estm_R=estm;
% % toc
% % break
%
%
% open ariport_final2.fig
% hold on
% plot(gd(:,2),gd(:,1),'o:')
% axis([126.45219 126.45255 37.44295 37.44325])
% xlabel('longitude')
% ylabel('latitude')
% title('161223 ��õ���� Pole 1')
% break

