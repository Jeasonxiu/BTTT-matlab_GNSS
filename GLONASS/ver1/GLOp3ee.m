%% GLONASS ��۱˵��� PPP �ڵ����� %% 2015.01.16
clear all; tic;
%% ��� ����
CCC = 299792458.;   % CCC = Speed of Light [m/s]
%% QM ����, ��¥ �ڵ鸵
DOY = 212; YY = 14;
% DOY = 225; YY = 16;
% --- GPS WEEK ����
[gw, gd] = ydoy2gwgd(YY, DOY);
% --- QM file �غ�
FileQM = 'IHUR2120';
% FileQM = 'uatA';
[arrQM, FinalPRNs, FinalTTs] = ReadQM(FileQM);

ObsType = 220; % ����ġ ����: GLONASS[211(L1), 212(L2), 220(C1), 221(P1), 222(P2), 231(D1), 232(D2)]
% ObsType = 320; % ����ġ ����: GLONASS[211(L1), 212(L2), 220(C1), 221(P1), 222(P2), 231(D1), 232(D2)]

QM = SelectQM(arrQM,ObsType);
FinalPRNs = unique(QM(:,2));
FinalTTs = unique(QM(:,1));
% --- ��ǥ �Է�
TruePos = [-3026676.0349  4067187.8095  3857246.8615 ];
% TruePos = [-3041235.578 4053941.677 3859881.013];   % : JPspace A point
%% ��� ���� �Է� �� ������ ����
% --- GLONASS ��۱˵��� ����: EPH, LeapSeceond, Tau_C
FileNavGLO = strcat('brdc', num2str(DOY), '0.', num2str(YY), 'g');
EphGlo = ReadEPH_glo(FileNavGLO);
Tau_c = ReadTauC(FileNavGLO);
% LeapSec = GetLeapSec(FileNavGLO);
% --- GPS ��۱˵��� ����: ������-Klobuchar ����
FileNav = strcat('brdc', num2str(DOY), '0.', num2str(YY), 'n');
[al, be] = GetALBE(FileNav);
LeapSec = GetLeapSec(FileNav);
% --- IONEX ����: DCB
FileIon = strcat('igsg', num2str(DOY), '0.', num2str(YY), 'i');
DCB = ReadDCB(FileIon);

%% ��ǥ �ʱ�ġ ���� �� ���浵 ��ȯ
AppPos = TruePos;
gd = xyz2gd(AppPos); AppLat = gd(1); AppLon = gd(2);
%% ������ ���� �Ű����� ����
Maxiter = 5;
EpsStop = 1e-4;
ctr = 1; eleCut=15;
x = [AppPos ctr]; x = x'; deltat = 15;
%%
NoEpochs = length(FinalTTs);
estm = zeros(NoEpochs, 5);
nEst = 0;

for j = 1:NoEpochs
    for iter = 1:Maxiter
        HTH = zeros(4,4);
        HTy = zeros(4,1);
        
        gs = FinalTTs(j);
        indexQM = find(QM(:,1) == gs);
        QM1e = QM(indexQM,:);
        NoSats = length(QM1e);
        
        vec_site = x(1:3)';
        % --- �������� ����� ������ ���; ����ũ�� �� �� ���
        ZHD = TropGPTh(vec_site, gw, gs); % GPT
        
        for i = 1:NoSats
            prn = QM1e(i,2);
            obs = QM1e(i,4);
            icol = PickEPH_GLO(EphGlo, prn, gs);
            
            TauN = EphGlo(icol,12); % : tau & gamma - �ð���� ������ �ʿ�
            GammaN = EphGlo(icol,13);
            ch_num = EphGlo(icol,16); % : channel number - ������ ������ �ʿ�
            
            %%
            STT = GetSTTbrdc_GLO(gs, EphGlo, icol, x(1:3), deltat); % ��ȣ���޽ð� ���
            tc = gs - STT; % ��ȣ���޽ð� ����
            tc = tc - LeapSec + Tau_c; % : LeapSecond & TauC ����
            
            % ��ȣ���޽ð� & LeapSecond ������ ���� ��ġ ����
            [sat_pos,sat_vel] = GetSatPosGLO(EphGlo,icol,tc,deltat);
            SatPos = sat_pos;
%             SatPos = PZ2WGS(sat_pos);
            SatPos = RotSatPos(SatPos,STT); % ���� ����ȿ�� ���
%             SatPos = SatPos';
            
            DistXYZ = SatPos - vec_site;
            DistNorm = norm(DistXYZ);
            
            [az,el] = xyz2azel(DistXYZ, AppLat, AppLon);
            %             if el>=eleCut
            %%
            dIono = Klo_R(vec_site, al, be, tc, SatPos,ch_num); % ������ ����
            dTrop = ZHD2SHD(gw, gs, vec_site, el, ZHD); % ����� ����
            
            SatVel = sat_vel;
            dRel = (-2/CCC^2) * dot(SatPos, SatVel); % ������ ȿ��
            dDCB = AppDCB_glo(DCB,prn);% DCB
            %%
            tsv = tc;
            tb = EphGlo(icol,2) + LeapSec; % GPStime - 16; / GLOtime + 16; Ȯ���ϱ�
            dtSat = TauN - GammaN*(tsv-tb) + Tau_c + dRel + dDCB;
            com = DistNorm + x(4) - CCC*dtSat + dTrop + dIono;
            %%
            y = obs - com;
            
            H = [-DistXYZ(1)/DistNorm -DistXYZ(2)/DistNorm -DistXYZ(3)/DistNorm 1];
            HTH = HTH + H'*H;
            HTy = HTy + H'*y;
            %             end  % if el>=eleCut
        end % --- for i = 1:NoSats
        xhat = inv(HTH) * HTy;
        x = x + xhat;
        if norm(xhat) < EpsStop;
            nEst = nEst + 1;
            estm(nEst,1) = gs;
            estm(nEst,2:5) = x(1:4);
            fprintf('%8d : %8.2f\n', j, norm(x(1:3)' - TruePos))
            break;
        end
    end % --- for iter = 1:Maxiter
end % % --- for j = 1:NoEpochs

% figure(200)
% plot(estm(:,1),estm(:,5),'-ob'); grid on;
toc; 

%% Analysis of Positioning Errors
estm = estm(1:nEst, :);
% [dXYZ, dNEV]=PosErrors(estm(:,1), TruePos,estm(:,2:4));
[dXYZ, dNEV] = PosTErrors2(estm(:,1), TruePos, estm(:,2:5));