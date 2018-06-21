
clc;
clear;

%% 변수 설정
CCC = 299792458;
ObsType = 120;
dtr = pi/180;
eleCut = 15;
YY=16;
DOY=75;
[gw, gd] = ydoy2gwgd(YY, DOY); 
% FileNav = 'brdc0750.16n';
FileNav = 'brdc0490.16n';
% FileQM = 'QMjprt075a.16o';
FileQM = 'QDBUU049_obr';
TruePos = [-3041241.741  4053944.143 3859873.640]; %: 대성디폴리스 B지점
% TruePos  = [-3026795.499 4067267.161 3857084.459]; %jprt
[arrQM, FinalPRNs, FinalTTs] = ReadQM(FileQM);
QM = SelectQM(arrQM, ObsType);
eph = ReadEPH(FileNav);
[al,be] = GetALBE(FileNav);                          
% AppPos  = App_pos('DBUU049r.16.obs');              
AppPos  = TruePos;
gd = xyz2gd(AppPos); AppLat = gd(1); AppLon = gd(2); 
%% 초기화
MaxIter = 4;
EpsStop = 1e-5;
ctr = 1;
x = [AppPos ctr]; x = x';
NoEpochs = length(FinalTTs);
nEst=0;
NoPara = 4;

% P 초기화 
P = eye(NoPara);
P(1:3,1:3) = P(1:3,1:3)*1;
P(4,4) = P(4,4)*0.5; 

% Q 초기화
Q = eye(NoPara);
Q(4,4) = 3000000;

for i=1:NoEpochs % 관측데이터 개수만큼 실행
    
     R = 0; indxHR = 0;    
     obs_a = []; com_a = [];  H_a=[];
%      H_a = zeros(1, NoPara); 
     indexQM = find(QM(:,1) == FinalTTs(i));
       QM_1 = QM(indexQM,:);
       NoSats = length(QM_1);
       gs = QM_1(1,1);
       vec_site = x(1:3)';
       ZHD = TropGPTh(vec_site, gw, gs);
       
              
     for sat=1:NoSats % 위성마다
            prn = QM_1(sat,2);
            obs = QM_1(sat,4);
            icol = PickEPH(eph,prn,gs);
            toe = eph(icol, 8); a = eph(icol, 19); b = eph(icol, 20); c = eph(icol, 21); Tgd = eph(icol, 23);
            STT = GetSTTbrdc(gs, prn, eph, vec_site);
            tc = gs - STT;
            vec_sat = GetSatPosNC(eph, icol, tc);                   
            vec_sat = RotSatPos(vec_sat, STT); 
            vec_rho = vec_sat - vec_site';
            rho = norm(vec_rho);  
            gd = xyz2gd(vec_site); AppLat = gd(1); AppLon = gd(2);  
            [az,el] = xyz2azel(vec_rho,AppLat,AppLon);    
            
            if el>eleCut

            I = ionoKlob(al, be, gs, az, el, vec_site);
            T = ZHD2SHD(gw,gs,vec_site,el,ZHD);
            dRel = GetRelBRDC(eph,icol,tc);
            dtSat = a+b*(tc-toe)+ c*(tc-toe)^2 -Tgd +dRel;
            com = rho- CCC*dtSat + T + I ;  
            indxHR = indxHR + 1;
            R(indxHR, indxHR) = 1/sin(el*dtr);         
            H = [ -vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho  1];
             obs_a = vertcat(obs_a, obs');
             com_a = vertcat(com_a, com');
             H_a = vertcat(H_a, H);
            end 
     end     
     xp = x;
     Po = P;
     Pp = P + Q;   
     K = Pp*H_a'*inv(H_a*Pp*H_a' + R);
     x = xp + K*(obs_a - com_a);
     P = Pp - K*H_a*Pp;             
     nEst = nEst + 1;
     estm(nEst,1) = gs;
     estm(nEst,2:5) = x(1:4);        
     fprintf('%8d : %8.2f\n',i, norm(x(1:3)' - TruePos))
end
[dXYZ, dNEV] = PosTErrors3(estm, TruePos);
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