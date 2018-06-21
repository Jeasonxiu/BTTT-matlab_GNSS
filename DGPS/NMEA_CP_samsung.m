%
%   modified by Joonseong Gim, aug 7, 2016
%

clear all; close all;
%% �Һ� ���� ����: ���� �ӵ�, ����ġ
CCC = 299792458.;   % CCC = Speed of Light [m/s]
%% �Ӱ���� ����
eleCut = 15;
%% GPS/GLO ����ġ ����
sys_w = [0.5, 0.5];
%% NMEA ���Ϸ� ���� NMEAQM ����
POINT = 'A';                    % logginf Point      
DOY = '050'; YY = '16';         % Day of Year    
TIME = 'i';                     % Logging Hour(9a,10b,11c,12d,13e,14f,15g,16h,17i,18j,19k,
                                %             20l,21m,22n,23o,24p,1q,2r,3s,4t,5u,6v,7w,8x)  
DEVICE = 'tab';

NMEAfile = strcat('jpr',POINT,DOY,TIME,'_',DEVICE,'.txt');
NMEAQM = writeNMEA3(NMEAfile);
% NMEAQM = writeNMEA2(NMEAfile);

%% load a PRC file
% PRCfile = 'JPRT160114.t1';      % tab
% PRCfile = 'JPRT160215.t1';      % tab
PRCfile = 'JPRT160219.t1';      % tab
% PRCfile = 'JPRT160201.t1';      % S6  
% PRCfile = 'JPRT160202.t1';      % S6
% PRCfile = 'JPRT160203.t1';      % S6
% PRCfile = 'JPRT160215.t1';      % S6
% PRCfile = 'JPRT160219.t1';      % S6
% PRCfile = 'JPRT160222.t1';      % S6
[GPSPRC, GLOPRC, PRC_Sorted] = PRCsort(PRCfile, NMEAQM); 
%% logging ������ ���� ������ ����
if POINT == 'A'
    TruePos = [-3041235.578 4053941.677 3859881.013];   % : JPspace A point
elseif POINT == 'B'
    TruePos = [-3041241.741 4053944.143 3859873.640];   % : JPspace B point
end

%% NMEAQM �ڵ鸵
FinalTTs = intersect(unique(NMEAQM(:,1)), unique(PRC_Sorted(:,1))); 

%% �׹��޽����� �о�鿩�� ��ķ� �����ϰ�, Klobuchar �� ����, gps week, gps day ����
navfile = strcat('brdc',DOY,'0','.',YY,'n');
eph = ReadEPH(navfile);
FileNavGLO = strcat('brdc',num2str(DOY),'0.',num2str(YY),'g');
EphGlo = ReadEPH_GLO(FileNavGLO);
TauC = ReadTauC(FileNavGLO); %%
LeapSec = GetLeapSec(navfile); %%
[al, be] = GetALBE(navfile);
[gw, GD] = ydoy2gwgd(str2num(YY), str2num(DOY)); %: GPS WEEK ����


%% ������ �ʿ��� �ʱ�ġ ����
MaxIter = 4;
EpsStop = 1e-5;
ctr = 1; deltat = 1;

%% �������� ����
NoEpochs = length(FinalTTs);
EstPos = zeros(NoEpochs,5);
nEst = 0;
j=1;
estm = zeros(NoEpochs,6);

% load('samsung4.mat');
x_prc = zeros(4,1);
sys_w = [0.9, 0.1];

for j = 1:NoEpochs
% for j = 1:804
    gs = FinalTTs(j);
    
    indexNMEAQM = find(NMEAQM(:,1) == FinalTTs(j));     % NMEA ��ü ��Ŀ��� ���� gs�� �� �� ����
    NMEAQM_1 = NMEAQM(indexNMEAQM,:);                   % NMEA ��ü ��Ŀ��� ����� gs ���� ������ ����
    indexPRC = find(PRC_Sorted(:,1) == FinalTTs(j));   % PRC ��ü ��Ŀ��� ���� gs�� �� �� ����
    PRC_1 = PRC_Sorted(indexPRC,:);                    % PRC ��ü ��Ŀ��� ����� gs ���� ������ ����
    
    NoSats = length(NMEAQM_1(:,5));                     % ���� gs�� ������
    
    ddmm = [NMEAQM_1(1,2) NMEAQM_1(1,3) NMEAQM_1(1,4)]; % NMEA ��ü ��Ŀ��� ����� gs ���� Lat, Long, Alt ������ ����
    gd = ddmm2gd(ddmm);                                 % NMEA ��ü ��Ŀ��� ����� gs ���� ������ ����
    
    vec_site(j,:) = gd2xyz(gd);                              % Receiver xyz
%     x = [vec_site(j,:) ctr ctr]; x = x';
    x = [vec_site(j,:) ctr]; x = x';
    ZHD = TropGPTh(vec_site(j,:), gw, gs);                 %: TROP: GPT
    
    GpsNMEAQM = find(NMEAQM_1(:,5) < 32); 
%     GpsSats = length(GpsQM);     % GPS C1
    GloNMEAQM = find(NMEAQM_1(:,5) > 32); 
%     GloSats = length(GloQM);     % GLO C1
%     visiSat(j,1) = gs; visiSat(j,2) = NoSats; visiSat(j,3) = GpsSats; visiSat(j,4) = GloSats;
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
    
    for Iter = 1:MaxIter
        HTH = zeros(4,4);                   % PP
        HTy = zeros(4,1);                   % PP
        HT_prc = zeros(4,1);
        
        for i = 1:NoSats
            prn = NMEAQM_1(i,5); 
            if prn < 32
                prc = PRC_1(find(PRC_1(:,2) == prn+100),3);               % DGPS PRC
                if ~isempty(prc)
                    icol = PickEPH(eph, prn, gs);
                    toe = eph(icol, 8); a = eph(icol, 19); b = eph(icol, 20); c = eph(icol, 21); Tgd = eph(icol, 23);
                    %----- ��ȣ���޽ð� ���
                    STT = GetSTTbrdc(gs, prn, eph, x(1:3)');                % GPS PP
                    tc = gs - STT;                                          % GPS PP

                    %----- �����˵� ���
                    vec_sat = GetSatPosNC(eph, icol, tc);                   % GPS PP
                    vec_sat = RotSatPos(vec_sat, STT);                      %: GPS PP �������� ���

                    %----- ���� RHO ���� ���
                    vec_rho = vec_sat - vec_site(j,:);                            % GPS PP
                    rho = norm(vec_rho);                                    % GPS
                    [az,el] = xyz2azel(vec_rho, gd(1), gd(2));   
                    
                    if el >= eleCut %15
                        W = 1;
%                         W = MakeW_elpr(el);
                        
                        dRel = GetRelBRDC(eph, icol, tc);                   % GPS PP
                        dtSat = a + b*(tc - toe) + c*(tc - toe)^2 - Tgd + dRel;                 % GPS PP
                        dIono = ionoKlob(al, be, gs, az, el, vec_site(j,:));
                        dTrop_G = ZHD2SHD(gw, gs, vec_site(j,:), el, ZHD);                   % GPT model

%                         y =  + dIono + dTrop_G;                                         % GPS PP
                        y =  - dIono -dTrop_G;                                         % GPS PP
%                         cp = + dIono + dTrop_G + prc;
                        cp = + dIono + dTrop_G + prc;
                        
                        H = [ -vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 1];                  % GPS PP
                        HTH = HTH + H'*W*sys_w(1)*H;                         % GPS PP
                        HTy = HTy + H'*W*sys_w(1)*y;                         % GPS PP
                        HT_prc = HT_prc + H'*W*sys_w(1)*cp;                                            % Correction-Projection
                    end
                end
            else
                
                % --- �������� ����� ������ ���: ����ũ�� �� �� ���
                prn = prn - 64;
                prc = PRC_1(find(PRC_1(:,2) == prn+300),3);
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
                    
                    DistXYZ = SatPos - vec_site(j,:);                         % GLO PP
                    DistNorm = norm(DistXYZ);                           % GLO PP
                    [az,el] = xyz2azel(vec_rho, gd(1), gd(2)); 
                    %%
                    if el>=eleCut
                        % ������ ���� % ����� ����
                        W = 1;
%                         W = MakeW_elpr(el);
                        
                        ttc = tc - TauC;                                    % GLO PP
                        dIono = Klo_R(vec_site(j,:),al,be,ttc,SatPos,ch_num);
                        dTrop = ZHD2SHD(gw,gs,vec_site(j,:),el,ZHD);
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
                      
%                         y =  + dIono + dTrop;
                        y =  - dIono ;
%                         cp = + dIono + dTrop + prc;
                        cp = + dIono + dTrop + prc;
                        
                        H = [-DistXYZ(1)/DistNorm -DistXYZ(2)/DistNorm -DistXYZ(3)/DistNorm 1];               % GLO PP
                        HTH = HTH + H'*W*sys_w(2)*H;                     % GLO PP
                        HTy = HTy + H'*W*sys_w(2)*y;                     % GLO PP
                        HT_prc = HT_prc + H'*W*sys_w(2)*cp;                                            % Correction-Projection
                    end
                end
            end
        end
        
        xhat_prc =  inv(HTH) * HTy - inv(HTH) * HT_prc;                               % Correction-Projection
%         xhat_prc =   - inv(HTH) * HT_prc;                               % Correction-Projection
        x_prc = x + xhat_prc;     
        
        if Iter == 4;
            nEst = nEst + 1;
            estm(nEst,1) = gs;
            estm_prc(nEst,1) = gs;
            estm(nEst,2:5) = x(1:4);
            estm_prc(nEst,2:5) = x_prc(1:4);
%             fprintf('gs: %6.0f     %2.1f    %2.1f \n',gs,norm(TruePos'-x(1:3)),norm(TruePos'-x_d(1:3)));
            break;
        end
    end
end
[dXYZ, dNEV] = PosTErrors22(estm(:, 1), TruePos, estm(:, 2:5));
[dXYZ, dNEV] = PosTErrorsCP(estm(:, 1), TruePos, estm(:, 2:5), estm_prc(:, 2:5));
[prcdXYZ, prcdNEV] = PosTErrors2(estm_prc(:, 1), TruePos, estm_prc(:, 2:5));

standalone = [mean(sqrt(dNEV(:,1).^2+dNEV(:,2).^2)), mean(dNEV(:,3)), std(sqrt(dNEV(:,1).^2+dNEV(:,2).^2)), std(dNEV(:,3)),rms(sqrt(dNEV(:,1).^2+dNEV(:,2).^2)), rms(dNEV(:,3))]
cpcpcp = [mean(sqrt(prcdNEV(:,1).^2+prcdNEV(:,2).^2)), mean(prcdNEV(:,3)), std(sqrt(prcdNEV(:,1).^2+prcdNEV(:,2).^2)), std(prcdNEV(:,3)),rms(sqrt(prcdNEV(:,1).^2+prcdNEV(:,2).^2)), rms(prcdNEV(:,3))]

for mm = 1: length(estm(:,1))
    smartgd(mm,1:3) = xyz2gd(estm(mm,2:4));
    cpgd(mm,1:3) = xyz2gd(estm_prc(mm,2:4));
end