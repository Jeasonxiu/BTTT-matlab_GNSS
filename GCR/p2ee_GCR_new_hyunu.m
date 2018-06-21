%% GPS & GLONASS & BDS ���� ���� �ڵ� (combined only)
tic
clear all; close all;
warning off;
%% ��� ����
CCC = 299792458.;   % CCC = Speed of Light [m/s]
%% QM ����, ��¥ �ڵ鸵
FileQM = 'B_Pole_a';
DOY = 355;
YY  = 16;
[arrQM, FinalPRNs, FinalTTs] = ReadQM(FileQM);
yyS = num2str(YY,'%02d');
doyS = num2str(DOY,'%03d');
[gw, gd] = ydoy2gwgd(YY, DOY); %: GPS WEEK ����
[bw, bd] = ydoy2bwbd(YY, DOY); %: BDS WEEK ����
%% ���� ��ǥ
% TruePos = [-2279828.8529 5004706.5404 3219777.4631]; %: JFNG site log 121029
% TruePos = [-1507972.56 6195614.06 148488.00]; %: SIN1 site log 130502
% TruePos = [-3026789.236 4067255.523 3857097.106];H=1.75; % inha
% TruePos = [-3026789.236 4067255.523 3857098.106]; % inha
% TruePos = [-3607665.4303 4147868.2248 3223717.2617];H=0; %: GMSD site log 120225
% TruePos = [-3026795.499 4067267.161 3857084.459]; % jprt
% TruePos = [-3027000.110018999   4066977.239737222   3857193.232413093];H=1.73; % PNT1
% TruePos = [-3026981.240   4066915.092   3857269.071]; % PNT7 (���׳� ���� 1.87m ���)
% TruePos = [-3113998.590 3920042.831 3938650.176];H=0; % �� PPS2
% TruePos = [-3290252.410 4018635.993 3690049.378];H=0; % ��� PPS3
TruePos = [-3041241.741 4053944.143 3859873.640]; % B
% TruePos = Apply_H(TruePos,H);
%% Navigation file load
FileNav = strcat('brdm', num2str(DOY,'%03d'), '0.', num2str(YY,'%02d'), 'p');   %: Navigation RINEX file
[eph, trashPRN, trashT]=ReadEPH_all(FileNav);
TauC = ReadTauC2(FileNav);
ephGLO = eph;
ephGLO=ephGLO(ephGLO(:,1)<400&ephGLO(:,1)>300,:);ephGLO(:,1)=ephGLO(:,1)-300;
% LeapSec = GetLeapSec(FileNav);
LeapSec = 17;
%% ���� �ڵ鸵
LeapSecBDS = 14;           %: BDS ����(2006.1.1)���� 2�� ����
%% brdm���Ͽ��� al/be�� ���� ������ brdc���Ͽ��� ������
al=[0 0 0 0]; be=[0 0 0 0]; %140
%% ����� ����ġ ����
g_ObsType = 120; % gps C1
c_ObsType = 220; % bds C1
r_ObsType = 320; % bds C1
%% QM �غ�
QM = SelectQM_gcr(arrQM, g_ObsType, c_ObsType, r_ObsType);
FinalPRNs = unique(QM(:,2));
FinalTTs = unique(QM(:,1));
%% ��ǥ �ʱ�ġ ���� �� ���浵 ��ȯ
AppPos = TruePos;
gd = xyz2gd(AppPos); AppLat = gd(1); AppLon = gd(2);
%% ������ ���� �Ű����� ����
Maxiter = 5;
EpsStop = 1e-4;
ctr = 1; ctr2 = 1; ctr3 = 1;
eleCut = 15; deltat = 5;
x = [AppPos ctr ctr2 ctr3]';
%%
NoEpochs = length(FinalTTs);
estm = zeros(NoEpochs, 7);
GDOP = zeros(NoEpochs, 1);
PDOP = zeros(NoEpochs, 1);
HDOP = zeros(NoEpochs, 1);
VDOP = zeros(NoEpochs, 1);
SatPosArr_before=0;icolArr_before=zeros(24,2);
nEst = 0;
for j = 1:NoEpochs % 60:65%
    icolArr=zeros(24,2);
    SatPosArr=zeros(24,11);
    for iter = 1:Maxiter
        gs = FinalTTs(j);
        HTH = zeros(6,6);
        HTy = zeros(6,1);
        indexQM = find(QM(:,1) == gs);
        QM1 = QM(indexQM,:);
        NoSats = length(QM1(:,1));
        vec_site = x(1:3)';  
        ZHD = TropGPTh(TruePos, gw, gs); %: TROP: GPT
        for i = 1:NoSats
            prn = QM1(i,2);
            obs = QM1(i,4);
            if prn > 100 && prn < 300
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
            elseif prn > 300 && prn < 400
                icol = PickEPH_GLO2(ephGLO, prn-300, gs);
                icolArr(i,1)=prn; icolArr(i,2)=icol;
                TauN = ephGLO(icol,12); % : tau & gamma - �ð���� ������ �ʿ�
                GammaN = ephGLO(icol,13);
                ch_num = ephGLO(icol,16); % : channel number - ������ ������ �ʿ�
                jcol = find(SatPosArr_before(:,1)==prn);
                if isempty(jcol); jcol=0;end
                jcol=jcol(1);
                %% ������ġ ��� ��Ʈ
                if icol == icolArr_before(icolArr_before(:,1)==prn,2) % icol ��ȭ ������ ���� epoch���� ���
                    %                 disp(1)
                    STT = GetSTTbrdc_GLO_ver3(gs, SatPosArr_before, jcol, x(1:3), deltat); % ���� epoch������ ��ȣ���޽ð� ���
                    tc = gs - STT;
                    tc = tc - LeapSec - TauC;
                    [SatPos,SatVel,SatLS] = GetSatPosGLO_ver3(SatPosArr_before,jcol,tc,deltat); % ���� epoch������ ������ġ ���
                else % icol ��ȭ��(������ �Ǵ� ���ο� ��۱˵���) ��۱˵������� ���
                    %                 disp(2)
                    STT = GetSTTbrdc_GLO(gs, ephGLO, icol, x(1:3), deltat); % ��۱˵������� ��ȣ���޽ð� ���
                    tc = gs - STT;
                    tc = tc - LeapSec - TauC;
                    [SatPos,SatVel,SatLS] = GetSatPosGLO(ephGLO,icol,tc,deltat); % ��۱˵������� ������ġ ���
                end
                
                SatPosArr(i,1)=prn; SatPosArr(i,2) = tc; SatPosArr(i,3:5) = SatPos;
                SatPosArr(i,6:8) = SatVel; SatPosArr(i,9:11) = SatLS;
            end
            SatPos = RotSatPos(SatPos, STT); % ��������ȿ�� ���
            vec_rho = SatPos - vec_site';
            rho = norm(vec_rho);
            [az,el] = xyz2azel(vec_rho, AppLat, AppLon);
%             dTrop = 0;
            dTrop = ZHD2SHD(gw, gs, TruePos, el, ZHD); % ����� ����
            
            if el >=eleCut % �Ӱ����
                %                 W = 1;
                W = MakeW_elpr(el);
                if prn > 100 && prn < 200
                    %             dIono = 0;
                    dIono = ionoKlob(al, be, gs, az, el, x(1:3)); % �̿��� ����(Klobuchar ��)
                    dRel = GetRelBRDC(eph, icol, tc); % ��뼺ȿ��
                    dtSat = a + b*(tc - toe) + c*(tc - toe)^2 - Tgd + dRel;
                    com = rho + x(4) - CCC*dtSat + dTrop + dIono; % gps ��갪
                    H = [-vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 1 0 0];
                elseif prn > 200 && prn < 300
                    %             dIono = 0;
                    dIono = ionoKlob(al, be, gs, az, el, x(1:3)); % �̿��� ����(Klobuchar ��)
                    %             dIono = ionoCIM(al, be, gs, az, el, x(1:3));   % �̿��� ����(Compass Ionospheric Model)
                    dRel = GetRelBRDC(eph, icol, tc); % ��뼺ȿ��
                    dtSat = a + b*(tc - toe) + c*(tc - toe)^2 - Tgd + dRel;
                    com = rho + x(5) - CCC*dtSat + dTrop + dIono; % bds ��갪
                    H = [-vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 0 1 0];
                elseif prn >300 && prn < 400
                    %             dIono = 0;
                    dIono = Klo_R(vec_site, al, be, tc, SatPos,ch_num); % ������ ����
                    dRel = (-2/CCC^2) * dot(SatPos, SatVel); % ������ ȿ��
                    dDCB = 0;
                    tsv = tc;
                    tb = ephGLO(icol,2) + LeapSec;
                    dtSat = TauN - GammaN*(tsv-tb) + TauC + dRel + dDCB;
                    com = rho + x(6) - CCC*dtSat + dTrop + dIono;
                    H = [-vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 0 0 1];
                end
                y = obs - com;
                HTH = HTH + H'*W*H;
                HTy = HTy + H'*W*y;
%                 fprintf('%8d : %3d : %d :::: obs: %f ::: com: %f ::: y: %f\n', j, iter, prn,obs,com,y);
            end
%             fprintf('%8d :: %3d :: PRN: %d :: SatPos: %8.3f %8.3f %8.3f :: icol: %d\n', j, iter, prn, SatPos,icol);
        end
        P = inv(HTH);
        xhat = inv(HTH) * HTy;
%         fprintf('%d\n', xhat');
%         fprintf('%8d : %3d : %f %f %f %f %f\n', j, iter, HTy');
        x = x + xhat;
        if norm(xhat) < EpsStop;
            nEst = nEst + 1;
            estm(nEst,1) = gs;
            estm(nEst,2:7) = x(1:6);
            GDOP(nEst) = sqrt(P(1,1)+P(2,2)+P(3,3)+P(4,4));
            PDOP(nEst) = sqrt(P(1,1)+P(2,2)+P(3,3));
            HDOP(nEst) = sqrt(P(1,1)+P(2,2));
            VDOP(nEst) = P(3,3);
            fprintf('%8d : %3d : %8.2f : %8.2f\n', j, iter, x(1)' - TruePos(1), x(2)' - TruePos(2));
            break;
        end
    end   
end
GDOP = GDOP(1:nEst,:);
PDOP = PDOP(1:nEst,:);
HDOP = HDOP(1:nEst,:);
VDOP = VDOP(1:nEst,:);
%% �������� �м� & �׷��� �ۼ�
estm = estm(1:nEst, :);
ttPlot=estm(:, 1);EstPosT=estm(:, 2:7);
% PosTErrors(ttPlot, TruePos, EstPosT);

%% �Է¹��� �ڷ��� ũ�� �� ������ �ش��ϴ� ���浵 ����
NoPos = length(EstPosT);
gd = xyz2gd(TruePos); TrueLat = gd(1); TrueLon = gd(2); 
%% ������ XYZ�� ���� XYZ�� ���̰� ����
dXYZ = zeros(NoPos,3);
for k = 1:NoPos
    dXYZ(k,:) = EstPosT(k,1:3) - TruePos;
end
%% dXYZ�� dNEV�� ��ȯ
dNEV = xyz2topo(dXYZ, TrueLat, TrueLon);

%% �� ���к� RMS ���
dN = dNEV(:,1); dE = dNEV(:,2); dV = dNEV(:,3);
dNE = sqrt(dN.^2 + dE.^2);        rmsH = sqrt(mean(dNE.^2));
                                  rmsV = sqrt(mean(dV.^2));
d3 = sqrt(dN.^2 + dE.^2 + dV.^2); rms3 = sqrt(mean(d3.^2));
dT_g = EstPosT(:,4);
dT_b = EstPosT(:,5);
dT_r = EstPosT(:,6);

fprintf('RMS Values \n H:%8.2f \n V:%8.2f \n 3D:%7.2f\n', rmsH, rmsV, rms3)
%% �׷��� �׸���
rXY = max([abs(min(dNEV(:,1))) abs(max(dNEV(:,1))) abs(min(dNEV(:,2))) abs(max(dNEV(:,2)))]);
rXY = ceil(rXY);

figure(1)
plot(dE, dN,'b.'); 
hold on
axis([-30 30 -30 30]);
grid on; 
xlabel('\Delta E (meters)'); ylabel('\Delta N (meters)')

tHour = mod(ttPlot, 86400);
tHour = tHour/3600;

figure(2)
subplot(3,2,1:2)
plot(tHour, dN, '.:'); 
axis([min(tHour) max(tHour) -20 20]); 
grid on;
ylabel('\Delta N (meters)')

subplot(3,2,3:4)
plot(tHour, dE, '.:'); 
axis([min(tHour) max(tHour) -20 20]); 
grid on; 
ylabel('\Delta E (meters)')

subplot(3,2,5:6)
plot(tHour, dV, '.:'); 
axis([min(tHour) max(tHour) -30 30]); 
grid on;
xlabel('Hours'); ylabel('\Delta U (meters)')
%%
% figure(1)
% set(figure(1), 'Position', [100, 100, 1200, 300]);
% subplot(1,9,1:3)
% plot(dE, dN,'ko','Markersize',10); 
% hold on
% axis([-5 5 -5 5]); 
% % axis([-10 10 -10 10]); 
% grid on; 
% xlabel('\Delta E (meters)'); ylabel('\Delta N (meters)')
% 
% subplot(1,9,5:9)
% tHour = mod(ttPlot, 86400);
% tHour = tHour/3600;
% plot(tHour, dV, 'k.:'); 
% axis([min(tHour) max(tHour) -10 10]); 
% % axis([min(tHour) max(tHour) -20 20]); 
% grid on;
% xlabel('Hours'); ylabel('\Delta V (meters)')
%%
figure(3)
subplot(3,2,1:2)
plot(tHour, dT_g, '.:'); axis([min(tHour) max(tHour) min(dT_g) max(dT_g)]); grid on; 
ylabel('GPS delta t_r (m)')
xlabel('Hours');

subplot(3,2,3:4)
plot(tHour, dT_b, '.:'); axis([min(tHour) max(tHour) min(dT_b) max(dT_b)]); grid on; 
ylabel('BDS delta t_r (m)')
xlabel('Hours');

subplot(3,2,5:6)
plot(tHour, dT_r, '.:'); axis([min(tHour) max(tHour) min(dT_r) max(dT_r)]); grid on; 
ylabel('GLO delta t_r (m)')
xlabel('Hours');
%%
% select_Sat='GC';
% PlotQM_NoSat(FileQM, select_Sat)
% legend('Total','GPS','BDS')
figure(5)
% hold on
plot(tHour,PDOP, 'm.:')
hold on
%%
% figure(4)
% plot(tHour, PDOP, 'k')
% hold on
% plot(tHour, HDOP, 'r')
% plot(tHour, VDOP, 'm')
% xlabel('Hours')
% ylabel('DOP')
% legend('PDOP','HDOP','VDOP')
toc