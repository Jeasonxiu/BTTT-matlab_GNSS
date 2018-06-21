close all; clear all;
tic
%% 불변 변수 설정: 빛의 속도, 관측치
CCC = 299792458.;   % CCC = Speed of Light [m/s]
ObsType = 120;      % 사용할 관측치 설정 - 120: C/A = C1
%% 임계고도각 설정
eleCut = 15;

%% PRC load
PRC = load('prc150523gps_sort.txt');
% PRCfile = 'prc150523new.txt';
% [PRC] = oldprc(PRCfile);
%% QMfile load
QMfile = 'Qihur15143';
% QMfile = 'QM15143_ihur';

%% QM 파일 읽어들여서 행렬로 저장하고, 사용할 관측치 추출
[arrQM, FinalPRNs, FinalTTs] = ReadQM(QMfile); 
QM = SelectQM(arrQM, ObsType);

%% 항법메시지를 읽어들여서 행렬로 저장하고, Klobuchar 모델 추출
eph = ReadEPH('brdc1430.15n');

%% 라이넥스 파일에서 대략적인 관측소 좌표를 뽑아내고 위경도로 변환
TruePos = [-3026675.978 4067187.900 3857246.933]; % IHUR (or IHU3) TRIMBLE NETR5
AppPos = TruePos;
gd = xyz2gd(AppPos); AppLat = gd(1); AppLon = gd(2); 

%% 추정에 필요한 초기치 설정
MaxIter = 4;
EpsStop = 1e-5;
ctr = 1;
x = [AppPos ctr]; x = x';           % PP
x_c = [AppPos ctr]; x_c = x_c';     % DGPS
x_cp = [AppPos ctr]; x_cp = x_cp';     % DGPS-CP

%% 추정과정 시작
FinalTTs = unique(PRC(:,1));
NoEpochs = length(FinalTTs);
nEst = 0;
nEst_c = 0;     % DGPS

for j = 1:NoEpochs
    for Iter = 1:MaxIter
        HTH = zeros(4,4);
        HTH_c = zeros(4,4);     % DGPS
        HTH_cp = zeros(4,4);    % DGPS-CP
        HTy = zeros(4,1);
        HTy_c = zeros(4,1);     % DGPS
        HTy_cp = zeros(4,1);    % DGPS-CP
        HT_cp = zeros(4,1);      % DGPS-CP
        HT_ppcp = zeros(4,1);      % pp-CP
                
        indexQM = find(QM(:,1) == FinalTTs(j));
        QM_1 = QM(indexQM,:);
        indexPRC = find(PRC(:,1) == FinalTTs(j));   % PRC 전체 행렬에서 현재 gs의 해 값 추출
        PRC_1 = PRC(indexPRC,:);                    % PRC 전체 행렬에서 추출된 gs 행의 데이터 추출
        NoSats = length(QM_1);
        gs = QM_1(1,1);
              
        vec_site = x(1:3)';             
        vec_site_c = x_c(1:3)';     % DGPS
        vec_site_cp = x_cp(1:3)';     % DGPS-CP
                
        for i = 1:NoSats
            prn = QM_1(i,2);
            obs = QM_1(i,4);   
            prc = PRC_1(find(PRC_1(:,2) == prn),3);   % DGPS
            icol = PickEPH(eph, prn, gs);
            toe = eph(icol, 8); a = eph(icol, 19); b = eph(icol, 20); c = eph(icol, 21); Tgd = eph(icol, 23);
            
            %----- 신호전달시간 계산
            STT = GetSTTbrdc(gs, prn, eph, vec_site);
            STT_c = GetSTTbrdc(gs, prn, eph, vec_site_c);    % DGPS
            STT_cp = GetSTTbrdc(gs, prn, eph, vec_site_cp);  % DGPS-CP
            tc = gs - STT;
            tc_c = gs - STT_c;                               % DGPS
            tc_cp = gs - STT_c;                              % DGPS-CP
                        
            %----- 위성궤도 계산
            vec_sat = GetSatPosNC(eph, icol, tc);          
            vec_sat_c = GetSatPosNC(eph, icol, tc_c);         % DGPS
            vec_sat_cp = GetSatPosNC(eph, icol, tc_cp);       % DGPS-CP
            vec_sat = RotSatPos(vec_sat, STT);                %: 지구자전 고려  
            vec_sat_c = RotSatPos(vec_sat_c, STT_c);          %: DGPS 지구자전 고려
            vec_sat_cp = RotSatPos(vec_sat_c, STT_cp);        %: DGPS-CP 지구자전 고려
                                    
            %----- 최종 RHO 벡터 계산
            vec_rho = (vec_sat - vec_site)';
            vec_rho_c = (vec_sat_c - vec_site_c)';            % DGPS
            vec_rho_cp = (vec_sat_c - vec_site_cp)';          % DGPS-CP
            rho = norm(vec_rho);
            rho_c = norm(vec_rho_c);                          % DGPS
            rho_cp = norm(vec_rho_cp);                        % DGPS-CP
            [az,el] = xyz2azel(vec_rho', AppLat, AppLon);        
                        
            if el >= eleCut %15
                W = 1;
                dRel = GetRelBRDC(eph, icol, tc);
                dRel_c = GetRelBRDC(eph, icol, tc_c);                                           % DGPS
                dRel_cp = GetRelBRDC(eph, icol, tc_cp);                                         % DGPS-CP
                dtSat = a + b*(tc - toe) + c*(tc - toe)^2 - Tgd + dRel;
                dtSat_c = a + b*(tc_c - toe) + c*(tc_c - toe)^2 - Tgd + dRel_c;                 % DGPS
                dtSat_cp = a + b*(tc_cp - toe) + c*(tc_cp - toe)^2 - Tgd + dRel_cp;             % DGPS-CP
                com = rho + x(4) - CCC * dtSat;
                com_c = rho_c + x_c(4) - CCC * dtSat_c - prc;                                   % DGPS
                com_cp = rho_cp + x_cp(4) - CCC * dtSat_cp;                                     % DGPS-CP
                cp = -prc;                                                                      % DGPS-CP
                y = obs - com;
                y_c = obs - com_c;                                                              % DGPS
                y_cp = obs - com_cp;                                                            % DGPS-CP
                H = [ -vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 1];
                H_c = [ -vec_rho_c(1)/rho_c -vec_rho_c(2)/rho_c -vec_rho_c(3)/rho_c 1];         % DGPS
                H_cp = [ -vec_rho_cp(1)/rho_cp -vec_rho_cp(2)/rho_cp -vec_rho_cp(3)/rho_cp 1];  % DGPS-CP
                HTH = HTH + H'*W*H;
                HTH_c = HTH_c + H_c'*W*H_c;                                             % DGPS
                HTH_cp = HTH_cp + H_cp'*W*H_cp;                                         % DGPS-CP
                HTy = HTy + H'*W*y;
                HTy_c = HTy_c + H_c'*W*y_c;                                             % DGPS
                HTy_cp = HTy_cp + H_cp'*W*y_cp;                                         % DGPS-CP
                HT_cp = HT_cp + H_cp'*W*cp;                                             % DGPS-CP_PRC
                HT_ppcp = HT_ppcp + H'*W*cp;                                            % pp-CP_PRC
            end
                        
        end
        xhat = inv(HTH) * HTy;
        xhat_c = inv(HTH_c) * HTy_c;                                                    % DGPS
        xhat_cp = inv(HTH_cp) * HTy_cp - inv(HTH_cp) * HT_cp;                           % DGPS-CP
        xhat_ppcp = inv(HTH) * HTy - inv(HTH) * HT_ppcp;                                % pp-CP
        x = x + xhat;
        x_c = x_c + xhat_c;                                                             % DGPS
        x_cp = x_cp + xhat_cp;                                                          % DGPS-CP
        x_ppcp = x + xhat_ppcp;                                                         % pp-CP
        
%         if norm(xhat) < EpsStop;
        if Iter == 4;
            nEst = nEst + 1;
            estm(nEst,1) = gs;
            estm_c(nEst,1) = gs;
            estm_cp(nEst,1) = gs;
            estm_ppcp(nEst,1) = gs;
            estm(nEst,2:5) = x(1:4);
            estm_c(nEst,2:5) = x_c(1:4);
            estm_cp(nEst,2:5) = x_cp(1:4);
            estm_ppcp(nEst,2:5) = x_ppcp(1:4);
            break;
        end

    end
    
    user_gd(j,:) = xyz2gd(estm(j,2:4)); % user의 xyz값 gd 로 변환
    user_xyz(j,:) = estm(j,2:4);        % user xyz값 행렬로 변환 
  
end

%% 측위오차 분석 & 그래프 작성
estm = estm(1:nEst, :);
estm_c = estm_c(1:nEst, :);
% [dXYZ, DGPSdXYZ, dNEV, DGPSdNEV] = PosTErrorsDGPS(estm(:, 1), TruePos, estm(:, 2:5), estm_c(:, 2:5));
[dXYZ, dNEV] = PosTErrors_result(estm(:, 1), TruePos, estm(:, 2:5));
[DGPSdXYZ, DGPSdNEV] = PosTErrors_result(estm_c(:, 1), TruePos, estm_c(:, 2:5));
[CPdXYZ, CPdNEV] = PosTErrors_result(estm_cp(:, 1), TruePos, estm_cp(:, 2:5));
[ppCPdXYZ, ppCPdNEV] = PosTErrors_result(estm_ppcp(:, 1), TruePos, estm_ppcp(:, 2:5));
fid_out = fopen('estm15143.txt', 'w');
fid_out_c = fopen('DGPSestm15143.txt', 'w');
fid_out_cp = fopen('DGPSestm15143_CP.txt', 'w');
fid_out_ppcp = fopen('DGPSestm15143_ppCP.txt', 'w');
for k = 1:nEst
    fprintf(fid_out, '%8d %13.3f %13.3f %13.3f %8.3f \n', estm(k,: ));
    fprintf(fid_out_c, '%8d %13.3f %13.3f %13.3f %8.3f \n', estm_c(k,: ));
    fprintf(fid_out_cp, '%8d %13.3f %13.3f %13.3f %8.3f \n', estm_cp(k,: ));
    fprintf(fid_out_ppcp, '%8d %13.3f %13.3f %13.3f %8.3f \n', estm_ppcp(k,: ));
end
fclose(fid_out);
fclose(fid_out_c);
fclose(fid_out_cp);
fclose(fid_out_ppcp);

format long
[dNE dU d3] = RMS(dNEV); rmsdNEV = [dNE dU d3]
[DGPSdNE DGPSdU DGPSd3] = RMS(DGPSdNEV); rmsDGPSdNEV = [DGPSdNE DGPSdU DGPSd3]
[CPdNE CPdU CPd3] = RMS(CPdNEV); rmsCPdNEV = [CPdNE CPdU CPd3]
[ppCPdNE ppCPdU ppCPd3] = RMS(ppCPdNEV); rmsppCPdNEV = [ppCPdNE ppCPdU ppCPd3]

figure(98)
hold on; grid on;
axis square
plot3(dNEV(:,2), dNEV(:,1), dNEV(:,3),'o'); 
plot3(DGPSdNEV(:,2), DGPSdNEV(:,1), DGPSdNEV(:,3),'r.'); 
plot3(CPdNEV(:,2), CPdNEV(:,1), CPdNEV(:,3),'g^'); 
plot3(ppCPdNEV(:,2), ppCPdNEV(:,1), ppCPdNEV(:,3),'m*'); 

figure(99)
hold on; grid on;
axis square
plot3(DGPSdNEV(:,2), DGPSdNEV(:,1), DGPSdNEV(:,3),'kd'...
    , CPdNEV(:,2), CPdNEV(:,1), CPdNEV(:,3),'g^'...
    , ppCPdNEV(:,2), ppCPdNEV(:,1), ppCPdNEV(:,3),'m*'); 
legend('DGPS', 'DGPS-CP', 'PP-CP')
xlabel({['  rms DGPS \Delta NE = ', num2str(rmsDGPSdNEV(1))...
    ,'   rms DGPS \Delta V = ', num2str(rmsDGPSdNEV(2))...
    ,'   rms DGPS \Delta 3D = ', num2str(rmsDGPSdNEV(3))],...
    ['rms D-CP \Delta NE = ', num2str(rmsCPdNEV(1))...
    ,'   rms D-CP \Delta V = ', num2str(rmsCPdNEV(2))...
    '     rms D-CP \Delta 3D = ', num2str(rmsCPdNEV(3))],...
    ['   rms PP-CP \Delta NE = ', num2str(rmsppCPdNEV(1))...
    ,'   rms PP-CP \Delta V = ', num2str(rmsppCPdNEV(2))...
    ,'   rms PP-CP \Delta 3D = ', num2str(rmsppCPdNEV(3))]});