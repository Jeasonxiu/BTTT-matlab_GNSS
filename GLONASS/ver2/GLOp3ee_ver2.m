%%% GLONASS ��۱˵��� PPP �ڵ����� %%%
% April 10th, 2015, Mi-So Kim, ���������ڵ忡�� �ӵ� ������

clc; clear all; tic;
%% ��� ����
CCC = 299792458.;   % CCC = Speed of Light [m/s]
%% QM ����, ��¥ �ڵ鸵
DOY = 212; YY = 14;
% --- GPS WEEK ����
[gw, gd] = ydoy2gwgd(YY,DOY); %%
% --- QM File �غ�
FileQM = 'IHUR2120';
[arrQM, FinalPRNs, FinalTTs] = ReadQM(FileQM );
ObsType = 220; % ����ġ ����: [GLONASS] 211(L1), 212(L2), 220(C1), 221(P1), 222(P2),231(D1), 232(D2)
QM = SelectQM(arrQM,ObsType);
FinalPRNs = unique(QM(:,2));
FinalTTs = unique(QM(:,1));
%% ��ǥ �ʱ�ġ ���� �� ���浵 ��ȯ
% --- ��ǥ �Է�
TruePos = [-3026676.0349  4067187.8095  3857246.8615]; %IHUR1980.14o
AppPos = TruePos;
% --- ��ǥ ���浵 ��ȯ
gd = xyz2gd(AppPos); AppLat = gd(1); AppLon = gd(2);
%% ��� ���� �Է� �� ������ ����
% --- GLONASS ��۱˵��� ���� : EphGLO, LeapSecond, TauC
FileNavGLO = strcat('brdc',num2str(DOY),'0.',num2str(YY),'g');
EphGlo = ReadEPH_GLO(FileNavGLO);
TauC = ReadTauC(FileNavGLO); %%
% --- GPS ��۱˵��� ����: ������-Klobuchar
FileNav = strcat('brdc',num2str(DOY),'0.',num2str(YY),'n');
[al,be] = GetALBE(FileNav); %%
LeapSec = GetLeapSec(FileNav); %%
% --- IONEX ����: DCB
FileIon = strcat('igsg',num2str(DOY),'0.',num2str(YY),'i');
DCB = ReadDCB(FileIon);
%% ������ ���� �Ű����� ����
Maxiter = 10; EpsStop = 1e-5;
ctr = 1; eleCut=15; deltat = 1;
x = [AppPos ctr]; x = x';

%%
NoEpochs = length(FinalTTs);
EstPos = zeros(NoEpochs,5);
nEst = 0;
for j = 1:NoEpochs
    indexQM = find(QM(:,1) == FinalTTs(j));
    QM_1 = QM(indexQM,:); NoSats = length(QM_1);
    gs = QM_1(1,1); 
    vec_site = x(1:3)';
    visiSat(j,2) = NoSats; visiSat(j,1) = gs;
    %% 1�� �������� ����, 30�п� ���� ��ġ ���(��۱˵��� 15��, 45��/Forward,Backward)
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
    % --- �������� ����� ������ ���: ����ũ�� �� �� ���
    ZHD = TropGPTh(vec_site, gw, gs); % GPT
    %% 
    for iter = 1:Maxiter
        HTH = zeros(4,4);
        HTy = zeros(4,1);
        for i = 1:NoSats
            prn = QM_1(i,2);
            obs = QM_1(i,4);
            
            tc = gs - 16;
            icol=PickEPH_GLO2(EphGlo, prn, tc);
            
            TauN=EphGlo(icol,12); GammaN=EphGlo(icol,13); %: tau & gamma �ð���� ������ ���
            ch_num=EphGlo(icol,16); %: channel number ������ ������ ���
            
            % ��ȣ���޽ð� ���
            STT = GetSTTbrdcGLO2(Sat_ar,gs,prn,x(1:3));           
            % LeapSecond & ��ȣ���� �ð��� ������ ���� ��ġ ����
            [SatPos, SatVel] = SatPosLS_STT(Sat_ar,gs,prn,LeapSec,STT,TauC);

            % ���� ����ȿ�� ���
            SatPos = RotSatPos(SatPos,STT);
%             SatPos = SatPos';
            
            DistXYZ = SatPos - x(1:3)';
            DistNorm = norm(DistXYZ);           
            [az,el] = xyz2azel(DistXYZ, AppLat, AppLon);
            %%
%             if el>=eleCut
                % ������ ���� % ����� ����
                ttc = tc - TauC;
                dIono = Klo_R(vec_site,al,be,ttc,SatPos,ch_num); 
                dTrop = ZHD2SHD(gw,gs,vec_site,el,ZHD);
                % ����� ȿ�� (�������ӵ� �̿�)
                dRel = (-2/CCC^2) * dot(SatPos, SatVel);
                % DCB ���
                dDCB = AppDCB_glo(DCB,prn);
                %%
                % �����ð���� tsv, tb
                tsv = tc; tb = EphGlo(icol,2) + 16; % GPStime - 16; / GLOtime + 16; Ȯ���ϱ�
%                 dtSat = TauN - GammaN*(tsv-tb) + TauC + dRel + dDCB;
                dtSat = TauN - GammaN*(tsv-tb) + TauC + dRel ;
                
                com = DistNorm + x(4) - CCC*dtSat + dTrop + dIono;
                y = obs - com;
                H = [-DistXYZ(1)/DistNorm -DistXYZ(2)/DistNorm -DistXYZ(3)/DistNorm 1];
                HTH = HTH + H'*H;
                HTy = HTy + H'*y;
%             end
        end
        xhat = inv(HTH) * HTy;
        x = x + xhat;
        if norm(xhat) < EpsStop;
            nEst = nEst + 1;
            EstPos(nEst,1) = gs;
            EstPos(nEst,2:5) = x(1:4);
            fprintf('gs: %6.0f     %2.5f \n',gs,norm(TruePos'-x(1:3)));    
            break;
        end
    end
end

%% Analysis of Positioning Errors
EstPos = EstPos(1:nEst, :);
[dXYZ, dNEV]=PosErrors1(EstPos(:,1), TruePos, EstPos(:,2:4),visiSat(:,2));
toc;
[dXYZ, dNEV]=PosTErrors3(EstPos(:,1), TruePos, EstPos(:,2:5),visiSat);
