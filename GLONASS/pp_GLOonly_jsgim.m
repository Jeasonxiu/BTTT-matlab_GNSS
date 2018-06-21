%% Clearing everything before each run
clear all; close all;
%% ���� ����: ���� �ӵ�, ����ġ
CCC = 299792458.;   %: CCC = Speed of Light [m/s]
ObsType = 320;      %: GLONASS C1
LeapSec = 18;
eleCut  = 5;
%% ��¥�� ����Ʈ 
YY = 17; DOY = 025; QMfile = 'QM170125_A'; TruePos = [-3041235.578 4053941.677 3859881.013];       % �뼺 A point
%% 
FileNav = strcat('brdm', num2str(DOY,'%03d'), '0.', num2str(YY,'%02d'), 'p');   %: Navigation RINEX file
%% QM�� eph_glo�� �����ص���
% [arrQM, dummy, dummy] = ReadQM(QMfile);
% [QM, FinalPRNs, FinalTTs] = SelectQMFinal(arrQM, ObsType);
% [eph, trashPRN, trashT] = ReadEPH_all(FileNav);
% ephGLO = eph;
% ephGLO = ephGLO(ephGLO(:,1) < 400 & ephGLO(:,1) > 300,:);
load 'ws_ppGLO';
%% 
[gw, gd] = ydoy2gwgd(YY, DOY); %: GPS WEEK ����
TauC = ReadTauC2(FileNav);
%% ��ǥ �ʱ�ġ ���� �� ���浵 ��ȯ
AppPos = TruePos;
gds = xyz2gd(AppPos); AppLat = gds(1); AppLon = gds(2);
%% ������ ���� �Ű����� ����
Maxiter = 5;
EpsStop = 1e-4;
cdtr = 1; 
deltat = 60;
x = [AppPos cdtr]';         

NoEpochs = length(FinalTTs);
estm = zeros(NoEpochs, 5);
%%  ���� epoch�� ���� epoch�� ���ؼ�, icol�� ��ȭ�� �ִ��� �Ǵ��ϱ� ���� arr���� ��Ʈ �Դϴ�.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SatPosArr_before = zeros(24,25); %:?? ���ٰ� �������� ���ϴ� ����?
icolArr_before=zeros(24,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nEst = 0;

for j = 1:NoEpochs
    
    gs = FinalTTs(j);
    indexQM = find(QM(:,1) == gs);
    QM1e = QM(indexQM,:);
    NoSats = length(QM1e);
    
    for iter = 1:Maxiter
        
        HTH = zeros(4,4);
        HTy = zeros(4,1);
        NoSatsUsed = 0;
        
        if NoSats <= 4
            break
        end
        vec_site = x(1:3)';
        ZHD = TropGPTh(TruePos, gw, gs); %: TROP: GPT
        
        for i = 1:NoSats
            
            prn = QM1e(i,2);
            obs = QM1e(i,4);
            
            icol = PickEPH_GLO2(ephGLO, prn+300, gs); %: PRN��ȣ ����
            icolArr(i,1) = prn+300; 
            icolArr(i,2) = icol;
            TauN = ephGLO(icol,12); %: tau & gamma - �ð���� ������ �ʿ�
            GammaN = ephGLO(icol,13);
            jcol = find(SatPosArr_before(:,1) == prn+300); %:PRN��ȣ ����
            if isempty(jcol); jcol = 0; end
            jcol = jcol(1);
            %% ������ġ ��� ��Ʈ
            if icol == icolArr_before(icolArr_before(:,1)==prn+300,2) % icol ��ȭ ������ ���� epoch���� ���
                STT = 0.075;   
                tc = gs - STT;
                tc = tc - LeapSec - TauC;
                [SatPos,SatVel,SatLS] = GetSatPosGLO_ver3(SatPosArr_before,jcol,tc,deltat); % ���� epoch������ ������ġ ���
                SatPos = RotSatPos(SatPos, STT); % ��������ȿ�� ���
            else % icol ��ȭ��(������ �Ǵ� ���ο� ��۱˵���) ��۱˵������� ���
                STT = 0.075;
                tc = gs - STT;
                tc = tc - LeapSec - TauC;
                [SatPos,SatVel,SatLS] = GetSatPosGLO(ephGLO,icol,tc,deltat); % ��۱˵������� ������ġ ���
                SatPos = RotSatPos(SatPos, STT); % ��������ȿ�� ���
            end
            icolArr_before(i,:) = icolArr(i,:);
            SatPosArr_before(i,:) = ephGLO(icol,:);
            
            vec_rho = SatPos - vec_site';
            rho = norm(vec_rho);
            [az,el] = xyz2azel(vec_rho, AppLat, AppLon);
            dIono = 0;
            dTrop = ZHD2SHD(gw, gs, TruePos, el, ZHD); % ����� ����
            dRel = (-2/CCC^2) * dot(SatPos, SatVel); % ������ ȿ��
            tb = ephGLO(icol,2) + LeapSec;
            dtSat = TauN - GammaN*(tc-tb) + TauC + dRel;
            
            if el >= eleCut % �Ӱ����
                com = rho + x(4) - CCC*dtSat + dTrop + dIono;
                H = [-vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 1];
                y = obs - com;
                HTH = HTH + H'*H;
                HTy = HTy + H'*y;
                NoSatsUsed = NoSatsUsed + 1;
            end
        end
        if NoSatsUsed >= 4
            xhat = inv(HTH) * HTy;
            x = x + xhat;
            if norm(xhat) < EpsStop;
                
                nEst = nEst + 1;
                estm(nEst,1) = gs;
                estm(nEst,2:5) = x(1:4);
                estm(nEst,6:7) = [NoSats, NoSatsUsed];      % ���� ���� ���� ��� ����
                estm(nEst,8:10) = xyz2gd(x(1:3));           % ���� ecef�� gd�� ��ȯ
                fprintf('%8d : %3d : %8.2f : %8.2f : %8.2f\n', j, iter, x(1)' - TruePos(1), x(2)' - TruePos(2), x(3)' - TruePos(3));
                break;
            end
        else
            break;
        end
    end
end
%% �������� �м� & �׷��� �ۼ�
estm = estm(1:nEst, :);
[dXYZ, dNEV] = PosTErrors2(estm(:,1), TruePos, estm(:,2:5),estm(:,6:7));

% %% ��������plot
% figure(99)
% hold on; grid on;
% plot(estm(:,10),estm(:,9),'bo')
% plot_google_map;