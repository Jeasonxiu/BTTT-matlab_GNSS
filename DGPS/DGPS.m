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

%% 추정과정 시작
FinalTTs = unique(PRC(:,1));
NoEpochs = length(FinalTTs);
nEst = 0;
nEst_c = 0;     % DGPS

for j = 1:NoEpochs
    for Iter = 1:MaxIter
        HTH = zeros(4,4);
        HTH_c = zeros(4,4);     % DGPS
        HTy = zeros(4,1);
        HTy_c = zeros(4,1);     % DGPS
                
        indexQM = find(QM(:,1) == FinalTTs(j));
        QM_1 = QM(indexQM,:);
        indexPRC = find(PRC(:,1) == FinalTTs(j));   % PRC 전체 행렬에서 현재 gs의 해 값 추출
        PRC_1 = PRC(indexPRC,:);                    % PRC 전체 행렬에서 추출된 gs 행의 데이터 추출
        NoSats = length(QM_1);
        gs = QM_1(1,1);
              
        vec_site = x(1:3)';             
        vec_site_c = x_c(1:3)';     % DGPS
                
        for i = 1:NoSats
            prn = QM_1(i,2);
            obs = QM_1(i,4);   
            prc = PRC_1(find(PRC_1(:,2) == prn),3);   % DGPS
            icol = PickEPH(eph, prn, gs);
            toe = eph(icol, 8); a = eph(icol, 19); b = eph(icol, 20); c = eph(icol, 21); Tgd = eph(icol, 23);
            
            %----- 신호전달시간 계산
            STT = GetSTTbrdc(gs, prn, eph, vec_site);
            STT_c = GetSTTbrdc(gs, prn, eph, vec_site_c);   % DGPS
            tc = gs - STT;
            tc_c = gs - STT_c;                              % DGPS
                        
            %----- 위성궤도 계산
            vec_sat = GetSatPosNC(eph, icol, tc);          
            vec_sat_c = GetSatPosNC(eph, icol, tc_c);       % DGPS
            vec_sat = RotSatPos(vec_sat, STT);              %: 지구자전 고려  
            vec_sat_c = RotSatPos(vec_sat_c, STT_c);        %: DGPS 지구자전 고려
                                    
            %----- 최종 RHO 벡터 계산
            vec_rho = (vec_sat - vec_site)';
            vec_rho_c = (vec_sat_c - vec_site_c)';          % DGPS
            rho = norm(vec_rho);
            rho_c = norm(vec_rho_c);                        % DGPS
            [az,el] = xyz2azel(vec_rho, AppLat, AppLon);        
                        
            if el >= eleCut %15
                W = 1;
                dRel = GetRelBRDC(eph, icol, tc);
                dRel_c = GetRelBRDC(eph, icol, tc_c);                                   % DGPS
                dtSat = a + b*(tc - toe) + c*(tc - toe)^2 - Tgd + dRel;
                dtSat_c = a + b*(tc_c - toe) + c*(tc_c - toe)^2 - Tgd + dRel_c;         % DGPS
                com = rho + x(4) - CCC * dtSat;
                com_c = rho_c + x_c(4) - CCC * dtSat_c - prc;                           % DGPS
                y = obs - com;
                y_c = obs - com_c;                                                      % DGPS
                H = [ -vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 1];
                H_c = [ -vec_rho_c(1)/rho_c -vec_rho_c(2)/rho_c -vec_rho_c(3)/rho_c 1]; % DGPS
                HTH = HTH + H'*W*H;
                HTH_c = HTH_c + H_c'*W*H_c;                                             % DGPS
                HTy = HTy + H'*W*y;
                HTy_c = HTy_c + H_c'*W*y_c;                                             % DGPS
            end
                        
        end
        xhat = inv(HTH) * HTy;
        xhat_c = inv(HTH_c) * HTy_c;                                                    % DGPS
        x = x + xhat;
        x_c = x_c + xhat_c;                                                             % DGPS
        
%         if norm(xhat) < EpsStop;
        if Iter == 4;
            nEst = nEst + 1;
            estm(nEst,1) = gs;
            estm_c(nEst,1) = gs;
            estm(nEst,2:5) = x(1:4);
            estm_c(nEst,2:5) = x_c(1:4);
            break;
        end

    end
    
    user_gd(j,:) = xyz2gd(estm(j,2:4)); % user의 xyz값 gd 로 변환
    user_xyz(j,:) = estm(j,2:4);        % user xyz값 행렬로 변환 
  
end

%% 측위오차 분석 & 그래프 작성
estm = estm(1:nEst, :);
estm_c = estm_c(1:nEst, :);
[dXYZ, DGPSdXYZ, dNEV, DGPSdNEV] = PosTErrorsDGPS(estm(:, 1), TruePos, estm(:, 2:5), estm_c(:, 2:5));
fid_out = fopen('estm15143.txt', 'w');
fid_out_c = fopen('DGPSestm15143.txt', 'w');
for k = 1:nEst
    fprintf(fid_out, '%8d %13.3f %13.3f %13.3f %8.3f \n', estm(k,: ));
    fprintf(fid_out, '%8d %13.3f %13.3f %13.3f %8.3f \n', estm_c(k,: ));
end
fclose(fid_out);
fclose(fid_out_c);

figure(99)
hold on; grid on;
axis square
plot3(dNEV(:,2), dNEV(:,1), dNEV(:,3),'o'); 
plot3(DGPSdNEV(:,2), DGPSdNEV(:,1), DGPSdNEV(:,3),'.'); 