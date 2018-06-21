%
%   modified by Joonseong Gim, aug 7, 2016
%

clear all; close all;

%% QMfile ȣ��
FileQM = 'QDRUA_ob153_1';

%% YY DOY �Է�
YY = 16;, DOY = 153;

%% load a PRC file
PRCfile = 'JPRT160601.t1';

%% �Һ� ���� ����: ���� �ӵ�, ����ġ
CCC = 299792458.;   % CCC = Speed of Light [m/s]
ObsType = 120;      % ����� ����ġ ���� - 120: C/A = C1

%% �Ӱ���� ����
eleCut = 15;

%% ������ǥ ����
TruePos = [-3041235.578 4053941.677 3859881.013];   % : JPspace A point
% TruePos = [-3041241.741 4053944.143 3859873.640];   % : JPspace B point

%% rinex file �Է�
obsfile = '160524-ubx1.obs';

%% ������� QMfile�� obsevation Rinex ���� ���� ��¥, �ð����� ����
rename = renameQMfile(obsfile);
[gw, GD] = ydoy2gwgd(YY, DOY); %: GPS WEEK ����

%% ��� ���� �Է� �� ������ ����
% --- GLONASS ��۱˵��� ���� : EphGLO, LeapSecond, TauC
FileNavGLO = strcat('brdc',num2str(DOY),'0.',num2str(YY),'g');
EphGlo = ReadEPH_GLO(FileNavGLO);
TauC = ReadTauC(FileNavGLO); %%
% --- GPS ��۱˵��� ����: ������-Klobuchar
FileNav = strcat('brdc',num2str(DOY),'0.',num2str(YY),'n');
LeapSec = GetLeapSec(FileNav); %%
%% GPS QM ���� �о�鿩�� ��ķ� �����ϰ�, ����� ����ġ ����
[arrQM, FinalPRNs, FinalTTs] = ReadQM(FileQM);
QMGPS = SelectQM(arrQM, ObsType);       % GPS C1
QMGLO = SelectQM(arrQM, 320);           % GLO C1
FinalQM = [QMGPS; QMGLO];               % GPS/GLO C1
FinalTTs = unique(FinalQM(:,1));

[GPSPRC, GLOPRC, PRC_sorted] = PRCsort(PRCfile, FinalQM);

%% �׹��޽����� �о�鿩�� ��ķ� �����ϰ�, Klobuchar �� ����
eph = ReadEPH(FileNav);
[al, be] = GetALBE(FileNav);
%% ���̳ؽ� ���Ͽ��� �뷫���� ������ ��ǥ�� �̾Ƴ��� ���浵�� ��ȯ
% AppPos = [-3041241.741 4053944.143 3859873.640];   % : JPspace B point
AppPos = GetAppPos(obsfile);
if AppPos(1) == 0
    AppPos = TruePos;
end

gd = xyz2gd(AppPos); AppLat = gd(1); AppLon = gd(2);

%% ������ �ʿ��� �ʱ�ġ ����
MaxIter = 10;
EpsStop = 1e-5;
ctr = 1; deltat = 1;

%% �������� ����
NoEpochs = length(FinalTTs);
EstPos = zeros(NoEpochs,5);
nEst = 0;
j=1;
estm = zeros(NoEpochs,6);
visiSat = zeros(NoEpochs,2);
GPS0_c = zeros(NoEpochs,32);
GLOo_c = zeros(NoEpochs,24);

% load('DGPSCPtest.mat')
% TruePos = [-3032234.51900000,4068599.11300000,3851377.46900000];

x = [AppPos ctr]; x = x';
x_cp = [AppPos ctr]; x_cp = x_cp';


% for j = 1:NoEpochs
for j = 1:908
    gs = FinalTTs(j);
    indexQM = find(FinalQM(:,1) == gs);         % gs epoch�� QM
    indexPRC = find(PRC_sorted(:,1) == gs);      % gs epoch�� PRC(GPS/GLO)
    QM1 = FinalQM(indexQM,:);           % GPS/GLO C1 1 Epoch
    PRC1 = PRC_sorted(indexPRC,:);
    NoSats = length(QM1(:,1));
    vec_site = x(1:3)';
    vec_site_cp = x_cp(1:3)';
    GpsQM = find(QM1(:,3) == 120); GpsSats = length(GpsQM);     % GPS C1
    GloQM = find(QM1(:,3) == 320); GloSats = length(GloQM);     % GLO C1
    visiSat(j,1) = gs; visiSat(j,2) = NoSats; visiSat(j,3) = GpsSats; visiSat(j,4) = GloSats;
    sT= mod(gs,86400)/3600;
    sTh = floor(sT); sTm = sT - sTh;
    
    if j == 1 % ���۽ð� ���� ��ġ array
        %         [Sat_ar] = GetSatPosGLO_new(EphGlo,gs,deltat);
        [Sat_ar] = GetSatPosGLO_my(EphGlo,gs,deltat);
        %         fprintf('���� �ð� ��ġ ���\n');
    elseif (j ~= 1 && sTm == 0) % ������ ��
        clear Sat_ar
        %         [Sat_ar] = GetSatPosGLO_new(EphGlo,gs,deltat);
        [Sat_ar] = GetSatPosGLO_my(EphGlo,gs,deltat);
        %         fprintf('%d�� ���� ����\n',sTh);
    elseif (j ~= 1 && mod(sTm,0.5) == 0) % 30�� �� ��
        clear Sat_ar
        %         [Sat_ar] = GetSatPosGLO_new(EphGlo,gs,deltat);
        [Sat_ar] = GetSatPosGLO_my(EphGlo,gs,deltat);
        %         fprintf('%d�� 30�� ����\n',sTh);
    end
    
    ZHD = TropGPTh(vec_site, gw, gs);
    for Iter = 1:MaxIter
        HTH = zeros(4,4);                   % PP
        HTy = zeros(4,1);                   % PP
        HTy_cp = zeros(4,1);                 % DGNSS
        
        for i = 1:NoSats
            prn = QM1(i,2); type = QM1(i,3); obs = QM1(i,4);
            if type == 120
%                 prc = PRC1(find(PRC1(:,2) == prn+100),3);               % DGPS PRC
%                 if ~isempty(prc)
%                     icol = PickEPH(eph, prn, gs);
%                     toe = eph(icol, 8); a = eph(icol, 19); b = eph(icol, 20); c = eph(icol, 21); Tgd = eph(icol, 23);
%                     %----- ��ȣ���޽ð� ���
%                     STT = GetSTTbrdc(gs, prn, eph, x(1:3)');                % GPS PP
%                     tc = gs - STT;                                          % GPS PP
%                     %----- �����˵� ���
%                     vec_sat = GetSatPosNC(eph, icol, tc);                   % GPS PP
%                     vec_sat = RotSatPos(vec_sat, STT);                      %: GPS PP �������� ���
%                     %----- ���� RHO ���� ���
%                     vec_rho = vec_sat - x(1:3)';                            % GPS PP
%                     rho = norm(vec_rho);                                    % GPS
%                     [az,el] = xyz2azel(vec_rho, AppLat, AppLon);
%                     
%                     if el >= eleCut %15
%                         W = 1;
% %                         W = MakeW_elpr(el);
%                         
%                         dRel = GetRelBRDC(eph, icol, tc);                   % GPS PP
%                         dtSat = a + b*(tc - toe) + c*(tc - toe)^2 - Tgd + dRel;                 % GPS PP
%                         dIono = ionoKlob(al, be, gs, az, el, vec_site);                 % Klobuchar
%                         dTrop_G = ZHD2SHD(gw, gs, vec_site, el, ZHD);                   % GPT model
%                         com = rho + x(4)  - CCC * dtSat + dIono + dTrop_G;              % GPS PP(Klobuchar, GPT Model)
%                         cp = - dIono - dTrop_G -prc;                            % DGNSS_CP
%                         
%                         y = obs - com;                                          % GPS PP
%                         
%                         H = [ -vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 1 0];                  % GPS PP
%                         HTH = HTH + H'*W*sys_w(1)*H;                         % GPS PP
%                         HTy = HTy + H'*W*sys_w(1)*y;                         % GPS PP
%                         HTy_cp = HTy_cp + H'*W*sys_w(1)*cp;                  % DGPS_CP
%                     end
%                 end
            else
                
                % --- �������� ����� ������ ���: ����ũ�� �� �� ���
                prc = PRC1(find(PRC1(:,2) == prn+300),3);
                if ~isempty(prc)
                    tc = gs - LeapSec;
                    icol=PickEPH_GLO2(EphGlo, prn, tc);
                    
                    TauN=EphGlo(icol,12); GammaN=EphGlo(icol,13); %: tau & gamma �ð���� ������ ���
                    ch_num=EphGlo(icol,16); %: channel number ������ ������ ���
                    
                    % ��ȣ���޽ð� ���
                    STT = GetSTTbrdcGLO2(Sat_ar,gs,prn,x(1:3));                 % GLO PP
                    % LeapSecond & ��ȣ���� �ð��� ������ ���� ��ġ ����
                    [SatPos, SatVel] = SatPosLS_STT(Sat_ar,gs,prn,LeapSec,STT,TauC);            % GLO PP
                                        
                    % ���� ����ȿ�� ���
                    SatPos = RotSatPos(SatPos,STT);                     % GLO PP
                    %             SatPos = SatPos';
                    
                    DistXYZ = SatPos - x(1:3)';                         % GLO PP
                    DistNorm = norm(DistXYZ);                           % GLO PP
                    [az,el] = xyz2azel(DistXYZ, AppLat, AppLon);
                    %%
                    if el>=eleCut
                        % ������ ���� % ����� ����
                        W = 1;
%                         W = MakeW_elpr(el);
                        
                        ttc = tc - TauC;                                    % GLO PP
                        dIono = Klo_R(vec_site,al,be,ttc,SatPos,ch_num);
                        dTrop = ZHD2SHD(gw,gs,vec_site,el,ZHD);
                        % ����� ȿ�� (�������ӵ� �̿�)
                        dRel = (-2/CCC^2) * dot(SatPos, SatVel);            % GLO PP
                        % DCB ���
                        %                     dDCB = AppDCB_glo(DCB,prn-200);
                        %%
                        % �����ð���� tsv, tb
                        tsv = tc;                                           % GLO PP
                        tb = EphGlo(icol,2) + LeapSec; % GPStime - 16; / GLOtime + 16; Ȯ���ϱ�
                        %                 dtSat = TauN - GammaN*(tsv-tb) + TauC + dRel + dDCB;
                        dtSat = TauN - GammaN*(tsv-tb) + TauC + dRel;           % GLO PP
                        
                        com = DistNorm + x(4) - CCC*dtSat + dTrop + dIono;      % GLO PP
                        cp = - dTrop - dIono - prc;                             % DGLO_CP 
                        
                        y = obs - com;                                          % GLO PP
                        
                        H = [-DistXYZ(1)/DistNorm -DistXYZ(2)/DistNorm -DistXYZ(3)/DistNorm 1];               % GLO PP
                        HTH = HTH + H'*W*H;                     % GLO PP
                        HTy = HTy + H'*W*y;                     % GLO PP
                        HTy_cp = HTy_cp + H'*W*cp;              % DGLO_cp
                    end
                end
            end
        end
        
        xhat = inv(HTH) * HTy;                              % PP
        xhat_cp = - inv(HTH) * HTy_cp;                        % DGNSS_CP
        x = x + xhat;                                       % PP
        x_cp = x + xhat_cp;                                 % DGNSS_CP
        
        if norm(xhat) < EpsStop;
            nEst = nEst + 1;
            estm(nEst,1) = gs;                              % PP
            estm_cp(nEst,1) = gs;                            % DGNSS_CP
            estm(nEst,2:5) = x(1:4);                        % PP
            estm_cp(nEst,2:5) = x_cp(1:4);                      % DGNSS_CP
            fprintf('gs: %6.0f     %2.1f    %2.1f \n',gs,norm(TruePos'-x(1:3)),norm(TruePos'-x_cp(1:3)));
            break;
        end
    end
end
estm = estm(1:nEst,:);
%% �������� �м� & �׷��� �ۼ�
[dXYZ, dNEV] = PosTErrors4(estm(:,1), TruePos, estm(:,2:5),visiSat);
[dXYZ_cp, dNEV_cp] = PosTErrors4(estm_cp(:,1), TruePos, estm_cp(:,2:5),visiSat);
