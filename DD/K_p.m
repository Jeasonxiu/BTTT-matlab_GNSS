
clc;
clear;

%% 불변 변수 설정 : 빛의 속도, 관측치
CCC = 299792458;
ObsType = 120;
dtr = pi/180;
eleCut = 15;
FileNav = 'brdc0750.16n';
% WriteObs('DBUU011c.16o');
% FileQM = 'QDBUU011c.16o';
FileQM = 'Qjprt075a.16o';


% TruePos = [-3041241.741  4053944.143 3859873.640]; %: 대성디폴리스 B지점
TruePos  = [-3026795.499 4067267.161 3857084.459]; %jprt
[arrQM, FinalPRNs, FinalTTs] = ReadQM(FileQM);
QM = SelectQM(arrQM, ObsType);
eph = ReadEPH(FileNav);
[al,be] = GetALBE(FileNav);                          % Rinex nav 에서 alpha, beta 로드   
% AppPos  = GetAppPos('jprt075a.16o');                          % Rinex obs 에서 apppos 로드
AppPos  = [-3041241.741 4053944.143 3859873.640];
% 
% al=[ 0.1118D-07 -0.7451D-08 -0.5960D-07  0.1192D-06 ]; % 대성디폴리스
% be=[ 0.1167D+06 -0.2294D+06 -0.1311D+06  0.1049D+07 ]; % 대성디폴리스
% AppPos = TruePos; % 참값 좌표로 변환
gd = xyz2gd(AppPos); AppLat = gd(1); AppLon = gd(2); 
%% 추정에 필요한 초기치 설정
MaxIter = 4;
EpsStop = 1e-5;
ctr = 1;%3.575261153706439e+006;
x = [AppPos ctr]; x = x';
NoEpochs = length(FinalTTs);
nEst=0;

NoPara = 3 + 1 ; % 좌표, 수신기 시계오차

% P 초기화 
P = eye(NoPara);
P(1:3,1:3) = P(1:3,1:3)*1;
P(4,4) = P(4,4)*0.5; 

% Q 초기화
Q = eye(NoPara);
 Q(1:3,1:3) = Q(1:3,1:3)*0.0001;
%  Q(1:3,1:3) = Q(1:3,1:3)*0.00001;
% Q(1:3,1:3) = Q(1:3,1:3)*8.3;
%  Q(4,4) = Q(4,4)*3e6;
%  Q(4,4) = Q(4,4)*3.575261153706439e+006;
% Q=zeros(4,4);

for i=1:NoEpochs % 관측데이터 개수만큼 실행
    
   
     %% R/H_a/obs_a/com_a 초기화 - 매 epoch마다 초기화 필요
     R = 0; indxHR = 0;    
     obs_a = []; com_a = [];  H_a=[];
%      H_a = zeros(1, NoPara); 
     indexQM = find(QM(:,1) == FinalTTs(i));
       QM_1 = QM(indexQM,:);
       NoSats = length(QM_1);
       gs = QM_1(1,1);
       vec_site = x(1:3)';
%        ZHD = TropGPTh(vec_site, gw, gs);                 %: TROP: GPT
     for sat=1:NoSats % 위성마다
            prn = QM_1(sat,2);
            obs = QM_1(sat,4);
            icol = PickEPH(eph,prn,gs);
            toe = eph(icol, 8); a = eph(icol, 19); b = eph(icol, 20); c = eph(icol, 21); Tgd = eph(icol, 23);
            STT = GetSTTbrdc(gs, icol, eph, vec_site);
            tc = gs - STT;
            vec_sat = GetSatPosNC(eph, icol, tc);                   % 위성 궤도 계산
            vec_sat = RotSatPos(vec_sat, STT); 
            vec_rho = vec_sat - vec_site;
            rho = norm(vec_rho);  
            gd = xyz2gd(vec_site); AppLat = gd(1); AppLon = gd(2);  % 위경도
            [az,el] = xyz2azel(vec_rho,AppLat,AppLon);    
            
            if el>eleCut
%             I = ionoKlob_my(al, be, gs, az, el, vec_site); 
            I = ionoKlob(al, be, gs, az, el, vec_site);
            T = Hopfield(el, 11, vec_site, 9999);;  
%             I=0; T=0;
            dRel = GetRelBRDC(eph,icol,tc);
            dtSat = a+b*(tc-toe)+ c*(tc-toe)^2 -Tgd +dRel;

            com = rho- CCC*dtSat + T + I ;  % 계산된 값
            indxHR = indxHR + 1;
            R(indxHR, indxHR) = 1/sin(el*dtr);% 2.25 ; % 노이즈 설정(measurement noise)
         
         
         
         %% 누적 파트
             H = [ -vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho  1];
             obs_a = vertcat(obs_a, obs');
             com_a = vertcat(com_a, com');
             H_a = vertcat(H_a, H);
            end 
     end
     
     %% (1) 추정값과 오차공분산 예측
     xp = x;
     Po = P;
     Pp = P + Q;                          %:시스템 노이즈가 공분산에 더해지는 과정
   
    %% (2) 칼만 이득 계산   
     K = Pp*H_a'*inv(H_a*Pp*H_a' + R);  %: R = Measurement Noise  
                                        %: 현재 관측치의 노이즈가 칼만게인에 더해지는 과정
    %% (3) 추정값 계산
     x = xp + K*(obs_a - com_a);         %: P, K, x가 업데이트 됨    
     
    %% (4) 오차 공분산 계산
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