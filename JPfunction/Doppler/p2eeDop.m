clc; clear all;
%% p2�� ��۱˵��� ��� �ڵ�-�ǻ�Ÿ� Point Positioing �˰�������, ee�� �� �������� ��ǥ ���� ����
% -- Modifications --
% 9/20/2014: P3PR_v1�� P3PR_BRDC�� �����ؼ� ����
%% �Һ� ���� ����: ���� �ӵ�, ����ġ
CCC = 299792458.;   % CCC = Speed of Light [m/s]
ObsType = 120;      % ����� ����ġ ���� - 120: C/A = C1
%% �Ӱ���� ����
eleCut = 15;
%% DOY & YY�� ������ ������ ���� ��Ī�� ���� �ο��ϰ��� ��
YY = 15;
DOY = 187;
% Global_PM(YY, DOY);             %: �۷ι� ���� ����
[gw, gd] = doy2gwgd(DOY, YY); %: GPS WEEK ����
FileNav = strcat('brdc',num2str(DOY,'%03d'),'0.',num2str(YY,'%02d'),'n'); %: �׹� RINEX ���� ����
% FileGIM = strcat('igsg',num2str(DOY),'0.',num2str(YY),'i'); %: GIM IONEX ���� ����
%% QM ���� �ڵ鸵 
% FileQM = 'QIHU1_11169a';
% FileQM = 'QIHU3_11169a';
% FileQM = 'QLAMT_12101';
% FileQM = 'QJINJ_13174';
% FileQM = 'QUBLX_14199';
% FileQM = 'QUBLX_14275';
% FileQM = 'QSBSA_14275p';
% FileQM = 'QIBU2_14328';
% FileQM = 'QIBU4_14328';
% FileQM = 'QIBU4_14331';
% FileQM = 'smo_14199'; %: QUBLX_14199 ������
% FileQM = 'smo_14275'; %: QUBLX_14275 ������
% FileQM = 'smo_sbs275'; %: QSBSA_14275p ������
% FileQM = 'QMfile';
% FileQM = 'QIVN1_15015';
% FileQM = 'QA1TB_15187';
%% ��Ÿ ����: ����Ʈ ��ǥ ���� & GPS Week
% TruePos = [-3026675.182 4067188.475 3857246.893]; %: IHU1 11169
% TruePos = [-3026675.978 4067187.900 3857246.933]; %: IHU3 11169
% TruePos = [-3037066.771 4049234.641 3867838.238]; %: UBLOX �Ѱ�������� 14199
% TruePos = [-3037068.684 4049234.342 3867837.343]; %: UBLOX �Ѱ�������� 14275
% TruePos = [-3003051.372 4059906.090 3883100.144]; %: SBSA ��ȭ 14275
% TruePos = [-3026450.541 4067350.305 3857210.992]; %: IBU2 ���ϴ� ���� UBLOX@PT2 14328/331
% TruePos = [-3026444.731 4067347.320 3857218.340]; %: IBU4 ���ϴ� ���� UBLOX@PT4 14328
% TruePos = [-3026789.236 4067255.523 3857098.106]; %: IVN1 ��ó�� ���� ������ 15014~
TruePos = [-3032239.745 4068575.036 3851397.953]; %: TB@�۵� 15187;
%% GIM �� ����� ���� �غ�
% [Lat, Lon, TEC] = ReadGIM(FileGIM);
%% ������������� GOA ������� �����ϱ� ���� ��� ����
% trop = ReadTROPgoa('TIHU3_13174');
%% QM ���� �о�鿩�� ��ķ� �����ϰ�, ����� ����ġ ����
FileQM = 'QA1TB_15187';
[arrQM, FinalPRNs, FinalTTs] = ReadQM(FileQM);      %* Read QMfile
% arrQM = load('QA1TB_15187.dat');
TTs = arrQM(:,1);
PRNs = arrQM(:,2);
FinalTTs = unique(TTs);
FinalPRNs = unique(PRNs);

QM = SelectQM(arrQM, ObsType);
%% �׹��޽����� �о�鿩�� ��ķ� �����ϰ�, Klobuchar �� ����
eph = ReadEPH(FileNav);
[al, be] = GetALBE(FileNav);
%% ���̳ؽ� ���Ͽ��� �뷫���� ������ ��ǥ�� �̾Ƴ��� ���浵�� ��ȯ
% AppPos = GetAppPos(FileObs); 
AppPos = TruePos;
gd = xyz2gd(AppPos); AppLat = gd(1); AppLon = gd(2); 
%% ������ �ʿ��� �ʱ�ġ ����
MaxIter = 4;
EpsStop = 1e-5;
ctr = 1;
x = [AppPos ctr]; x = x';
%% �������� ����
NoEpochs = length(FinalTTs);
estm = zeros(NoEpochs, 8); % c6/7/8(Vx, Vy, Vz)
nEst = 0;

for kE = 481:NoEpochs
    indexQM = find(QM(:,1) == FinalTTs(kE));
    QM1e = QM(indexQM,:);
    NoSats1e = length(QM1e);
    
    gs = QM1e(1,1);
    

    
    nDop = 0;                       %** 7/12/15�߰� DOP �ӵ�
    prnsDop = zeros(NoSats1e, 1);
    dops = zeros(NoSats1e, 1);
    

    
    for Iter = 1:MaxIter
        HTH = zeros(4,4);
        HTy = zeros(4,1);
        
        vec_site = x(1:3)';
        ZHD = TropGPTh(vec_site, gw, gs);                 %: TROP: GPT
                
        for kS = 1:NoSats1e
            prn = QM1e(kS,2);
            obs = QM1e(kS,4);
            
            icol = PickEPH(eph, prn, gs);
            toe = eph(icol, 8); a = eph(icol, 19); b = eph(icol, 20); c = eph(icol, 21); Tgd = eph(icol, 23);
            
            %----- ��ȣ���޽ð� ���
            STT = GetSTTbrdc(gs, eph, icol, vec_site);
            tc = gs - STT;
            %----- �����˵� ���
            vec_sat = GetSatPosNC(eph, icol, tc);
            vec_sat = RotSatPos(vec_sat, STT);                  %: �������� ���
            %----- ���� RHO ���� ���
            vec_sat = vec_sat';
            vec_rho = vec_sat - vec_site;
            rho = norm(vec_rho);
            
            [az,el] = xyz2azel(vec_rho, AppLat, AppLon);        %: ���߿� ������ ��. �ʱ�ġ ���� �ʿ�
            
            if el >= eleCut
                
                if Iter == 1
                    nDop = nDop + 1;
                    indxDop = find(arrQM(:,1) == gs & arrQM(:,2) == prn & arrQM(:,3) == 131);
                    
                    prnsDop(nDop, 1) = prn;
                    dops(nDop, 1) = arrQM(indxDop, 4);
                end
                
                W = 1;
                % W = MakeW_elpr(el);
                dIono = ionoKlob(al, be, gs, az, el, vec_site);           %: IONO: Klobuchar
                % dIono = IonoGIM(Lat, Lon, TEC, eph, x(1:3)', prn, gs);    %: IONO: GIM
                dTrop = ZHD2SHD(gw, gs, vec_site, el, ZHD);
                
                dRel = GetRelBRDC(eph, icol, tc);
                dtSat = a + b*(tc - toe) + c*(tc - toe)^2 - Tgd + dRel;
                
                com = rho + x(4) - CCC * dtSat + dIono + dTrop;
                y = obs - com;
                H = [ -vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 1];
                HTH = HTH + H'*W*H;
                HTy = HTy + H'*W*y;
                
            end
        end
        xhat = inv(HTH) * HTy;
        x = x + xhat;
        
        if norm(xhat) < EpsStop;
            nEst = nEst + 1;
            estm(nEst,1) = gs;
            estm(nEst,2:5) = x(1:4);

            estm(nEst,6:8) = GetVelDop(gs, prnsDop(1:nDop), dops(1:nDop), eph, x(1:3)');
            % fprintf('%8d : %8.2f\n',j, norm(x(1:3)' - TruePos))
            break;
        end
        
    end
end
%% �������� �м� & �׷��� �ۼ�
estm = estm(1:nEst, :);
% [dXYZ, dNEV] = PosTErrors(estm(:, 1), TruePos, estm(:, 2:5));
[dXYZ, dNEV] = PosErrors(estm(:, 1), TruePos, estm(:, 2:4));
fid_out = fopen('estm.txt', 'w');
for k = 1:nEst
    fprintf(fid_out, '%8d %13.3f %13.3f %13.3f %8.3f \n', estm(k,: ));
end
fclose(fid_out);
%% Doppler �ӵ� �׷��� �ۼ�
figure(777)
PlotVelDop(estm);

