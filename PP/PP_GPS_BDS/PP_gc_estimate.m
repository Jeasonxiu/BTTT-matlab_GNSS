
clear all;
% %
% QMfile = 'QM170125_A';
% QMfile = 'T4_170323_k500';
QMfile = 'T4_GC';
% QMfile = 'QMfile';
% navfile = 'brdm0250.17p';
% navfile = 'hour0450.17n';
navfile = 'brdm0820.17p';
% TruePos = [-3041235.578 4053941.677 3859881.013];
TruePos = [-3052725.48122377,4042747.35562486,3862428.65171442];
YY = 17; DOY =082;
% TruePos = [-3058799.61420451,4083265.35912516,3814946.87192938];

%% �Һ� ���� ����: ���� �ӵ�, ����ġ
CCC = 299792458.;   % CCC = Speed of Light [m/s]
%% �Ӱ���� ����
eleCut = 15;
%% ���� �ڵ鸵
LeapSecBDS = 14;

%% ����� ����ġ ����
g_ObsType = 120; % gps C1
g_ObsType_snr = 141;
c_ObsType = 220; % bds C1
c_ObsType_snr = 241;

%% QM ���� �о�鿩�� ��ķ� �����ϰ�, ����� ����ġ ����
[arrQM, FinalPRNs, FinalTTs] = ReadQM(QMfile);
arrQM = arrQM(find(arrQM(:,3) < 300),:);
QM = SelectQM_gc(arrQM, g_ObsType, c_ObsType);
QM_snr = SelectQM_gc(arrQM, g_ObsType_snr, c_ObsType_snr);
FinalTTs = unique(QM(:,1));
% FinalTTs = FinalTTs(find(FinalTTs(:,1) > 1000),:);

% %% eph ����
[eph, trashPRN, trashT]=ReadEPH_all(navfile);

%% �׹��޽����� �о�鿩�� ��ķ� �����ϰ�, Klobuchar �� ����
gps_nav = strcat('brdm', num2str(DOY,'%03d'), '0.', num2str(YY,'%02d'), 'p');   %: Navigation RINEX file
fid = fopen(gps_nav,'r');
if fid == -1
    al = zeros(4,1); be = zeros(4,1);
else
    [al, be] = GetALBE(gps_nav);
end

%% �׹��޽��� ���� �̸��� �̿��� YY, DOY ����
% YY = str2num(navfile(length(navfile)-2:length(navfile)-1));
% DOY = str2num(navfile(length(navfile)-7:length(navfile)-5));
[gw, gd] = ydoy2gwgd(YY, DOY); %: GPS WEEK ����

%% �ʱ���ǥ ȹ��
AppPos = TruePos;
gd = xyz2gd(AppPos); AppLat = gd(1); AppLon = gd(2);

%% ������ �ʿ��� �ʱ�ġ ����
Maxiter = 10;
EpsStop = 1e-5;
ctr = 1; ctr2 = 1;
x = [AppPos ctr ctr2]; x = x';

%% �������� ����
NoEpochs = length(FinalTTs);
estm = zeros(NoEpochs, 6);
Scount = zeros(NoEpochs, 1);
nEst = 0;

for j = 1:NoEpochs
    gs = FinalTTs(j);
    indexQM = find(QM(:,1) == gs);
    QM1e = QM(indexQM,:);
    QM1e_snr = QM_snr(indexQM,:);
    NoSats = length(QM1e);
   
    for iter = 1:Maxiter
        HTH = zeros(5,5);
        HTy = zeros(5,1);
        
        if NoSats <= 6
            break
        end
        vec_site = x(1:3)';
        ZHD = TropGPTh(TruePos, gw, gs); %: TROP: GPT
        for i = 1:NoSats
            NoSatsUsed = 0;
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
                %                     dIono = 0;
                %                     dTrop = 0;
                dIono = ionoKlob(al, be, gs, az, el, x(1:3)); % �̿��� ����(Klobuchar ��)
                dTrop = ZHD2SHD(gw, gs, TruePos, el, ZHD); % ����� ����
            elseif prn > 200 && prn < 300
                %                     dIono = 0;
                %                     dTrop = 0;
                dIono = ionoKlob(al, be, gs, az, el, x(1:3)); % �̿��� ����(Klobuchar ��)
                dTrop = ZHD2SHD(gw, gs, TruePos, el, ZHD); % ����� ����
            end
            %             dRel=0;
            dRel = GetRelBRDC(eph, icol, tc); % ��뼺ȿ��
            dtSat = a + b*(tc - toe) + c*(tc - toe)^2 - Tgd + dRel; % �����ð���� ���... �׷������, ��뼺ȿ�� ����
            if el >=eleCut % �Ӱ����
                %                 W = 1;
                W = MakeW_elsnr(el,snr);
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
            Scount(nEst,1) = NoSatsUsed;
            fprintf('%8d : %3d : %8.2f : %8.2f\n', j, iter, x(1)' - TruePos(1), x(2)' - TruePos(2));
            break;
        end
    end
end
estm = estm(1:nEst,:);
Scount=Scount(1:nEst, :);
