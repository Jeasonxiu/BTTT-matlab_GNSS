
clc;
clear;

%% �Һ� ���� ���� : ���� �ӵ�, ����ġ
CCC = 299792458;
ObsType = 120;
dtr = pi/180;
eleCut = 15;
FileNav = 'brdc0750.16n';
% WriteObs('DBUU011c.16o');
% FileQM = 'QDBUU011c.16o';
FileQM = 'Qjprt075a.16o';


% TruePos = [-3041241.741  4053944.143 3859873.640]; %: �뼺�������� B����
TruePos  = [-3026795.499 4067267.161 3857084.459]; %jprt
[arrQM, FinalPRNs, FinalTTs] = ReadQM(FileQM);
QM = SelectQM(arrQM, ObsType);
eph = ReadEPH(FileNav);
[al,be] = GetALBE(FileNav);                          % Rinex nav ���� alpha, beta �ε�   
% AppPos  = GetAppPos('jprt075a.16o');                          % Rinex obs ���� apppos �ε�
AppPos  = [-3041241.741 4053944.143 3859873.640];
% 
% al=[ 0.1118D-07 -0.7451D-08 -0.5960D-07  0.1192D-06 ]; % �뼺��������
% be=[ 0.1167D+06 -0.2294D+06 -0.1311D+06  0.1049D+07 ]; % �뼺��������
% AppPos = TruePos; % ���� ��ǥ�� ��ȯ
gd = xyz2gd(AppPos); AppLat = gd(1); AppLon = gd(2); 
%% ������ �ʿ��� �ʱ�ġ ����
MaxIter = 4;
EpsStop = 1e-5;
ctr = 1;%3.575261153706439e+006;
x = [AppPos ctr]; x = x';
NoEpochs = length(FinalTTs);
nEst=0;

NoPara = 3 + 1 ; % ��ǥ, ���ű� �ð����

% P �ʱ�ȭ 
P = eye(NoPara);
P(1:3,1:3) = P(1:3,1:3)*1;
P(4,4) = P(4,4)*0.5; 

% Q �ʱ�ȭ
Q = eye(NoPara);
 Q(1:3,1:3) = Q(1:3,1:3)*0.0001;
%  Q(1:3,1:3) = Q(1:3,1:3)*0.00001;
% Q(1:3,1:3) = Q(1:3,1:3)*8.3;
%  Q(4,4) = Q(4,4)*3e6;
%  Q(4,4) = Q(4,4)*3.575261153706439e+006;
% Q=zeros(4,4);

for i=1:NoEpochs % ���������� ������ŭ ����
    
   
     %% R/H_a/obs_a/com_a �ʱ�ȭ - �� epoch���� �ʱ�ȭ �ʿ�
     R = 0; indxHR = 0;    
     obs_a = []; com_a = [];  H_a=[];
%      H_a = zeros(1, NoPara); 
     indexQM = find(QM(:,1) == FinalTTs(i));
       QM_1 = QM(indexQM,:);
       NoSats = length(QM_1);
       gs = QM_1(1,1);
       vec_site = x(1:3)';
%        ZHD = TropGPTh(vec_site, gw, gs);                 %: TROP: GPT
     for sat=1:NoSats % ��������
            prn = QM_1(sat,2);
            obs = QM_1(sat,4);
            icol = PickEPH(eph,prn,gs);
            toe = eph(icol, 8); a = eph(icol, 19); b = eph(icol, 20); c = eph(icol, 21); Tgd = eph(icol, 23);
            STT = GetSTTbrdc(gs, icol, eph, vec_site);
            tc = gs - STT;
            vec_sat = GetSatPosNC(eph, icol, tc);                   % ���� �˵� ���
            vec_sat = RotSatPos(vec_sat, STT); 
            vec_rho = vec_sat - vec_site;
            rho = norm(vec_rho);  
            gd = xyz2gd(vec_site); AppLat = gd(1); AppLon = gd(2);  % ���浵
            [az,el] = xyz2azel(vec_rho,AppLat,AppLon);    
            
            if el>eleCut
%             I = ionoKlob_my(al, be, gs, az, el, vec_site); 
            I = ionoKlob(al, be, gs, az, el, vec_site);
            T = Hopfield(el, 11, vec_site, 9999);;  
%             I=0; T=0;
            dRel = GetRelBRDC(eph,icol,tc);
            dtSat = a+b*(tc-toe)+ c*(tc-toe)^2 -Tgd +dRel;

            com = rho- CCC*dtSat + T + I ;  % ���� ��
            indxHR = indxHR + 1;
            R(indxHR, indxHR) = 1/sin(el*dtr);% 2.25 ; % ������ ����(measurement noise)
         
         
         
         %% ���� ��Ʈ
             H = [ -vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho  1];
             obs_a = vertcat(obs_a, obs');
             com_a = vertcat(com_a, com');
             H_a = vertcat(H_a, H);
            end 
     end
     
     %% (1) �������� �������л� ����
     xp = x;
     Po = P;
     Pp = P + Q;                          %:�ý��� ����� ���л꿡 �������� ����
   
    %% (2) Į�� �̵� ���   
     K = Pp*H_a'*inv(H_a*Pp*H_a' + R);  %: R = Measurement Noise  
                                        %: ���� ����ġ�� ����� Į�����ο� �������� ����
    %% (3) ������ ���
     x = xp + K*(obs_a - com_a);         %: P, K, x�� ������Ʈ ��    
     
    %% (4) ���� ���л� ���
     P = Pp - K*H_a*Pp;       

            nEst = nEst + 1;
            estm(nEst,1) = gs;
            estm(nEst,2:5) = x(1:4);        
            fprintf('%8d : %8.2f\n',i, norm(x(1:3)' - TruePos))
end
[dXYZ, dNEV, dNE, NEV] = PosTErrors3(estm, TruePos);
figure(1)
subplot(1,2,1)
plot(dXYZ(:,1),dXYZ(:,2),'o')
lineCross(0,0,'r',2);
grid on
subplot(3,2,2)
plot(dNE(:,1)); legend('2d');
grid on
subplot(3,2,4)
plot(dXYZ(:,3)); legend('h');
grid on
subplot(3,2,6)
plot(NEV(:,1)); legend('3d');
grid on
dNEV