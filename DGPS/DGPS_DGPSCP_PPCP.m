close all; clear all;
tic
%% 불변 변수 설정: 빛의 속도, 관측치
CCC = 299792458.;   % CCC = Speed of Light [m/s]
ObsType = 120;      % 사용할 관측치 설정 - 120: C/A = C1
%% 임계고도각 설정
eleCut = 15;

%% gps week, gps day, alpha, beta 생성
% YY = 16; DOY = 049;
YY = 16; DOY = 145;
[gw, GD] = ydoy2gwgd(YY, DOY); %: GPS WEEK 결정
% navfile = 'brdc0490.16n';
navfile = 'brdc1450.16n';
[al, be] = GetALBE(navfile);
%% QMfile load
% QMfile = 'Qihur15143';
% QMfile = 'QMDBUU049i_16';
QMfile = 'Q160524_ubx1';
% QM = load(QMfile);
QM = load('Q160524_ubx1');
% QMfile = 'QM15143_ihur';
load('160524_1_adm.txt');
%% QM 파일 읽어들여서 행렬로 저장하고, 사용할 관측치 추출
[arrQM, FinalPRNs, FinalTTs] = ReadQM(QMfile); 
QM = SelectQM(arrQM, ObsType);

%% PRC load
% PRCfile = 'JPRT160218.t1';
PRCfile = 'JPRT160524.t1';
[PRC_Sorted] = PRCsort(PRCfile, QM); 

%% 항법메시지를 읽어들여서 행렬로 저장하고, Klobuchar 모델 추출
% eph = ReadEPH('brdc0490.16n');
eph = ReadEPH(navfile);

%% 라이넥스 파일에서 대략적인 관측소 좌표를 뽑아내고 위경도로 변환
% TruePos = [-3026675.978 4067187.900 3857246.933]; % IHUR (or IHU3) TRIMBLE NETR5
% TruePos = [-3041235.578 4053941.677 3859881.013];   % : JPspace A point
% TruePos = [-3041241.741 4053944.143 3859873.640];   % : JPspace B point
% AppPos = TruePos;

AppPos = GetAppPos('160524-ubx1.obs');
gd = xyz2gd(AppPos); AppLat = gd(1); AppLon = gd(2); 

%% 추정에 필요한 초기치 설정
MaxIter = 4;
EpsStop = 1e-5;
ctr = 1;
x = [AppPos ctr]; x = x';           % PP
x_c = [AppPos ctr]; x_c = x_c';     % DGPS
x_ppcp = [AppPos ctr]; x_ppcp = x_ppcp';
%% 추정과정 시작
% FinalTTs = unique(PRC(:,1));
FinalTTs = intersect(unique(QM(:,1)), unique(PRC_Sorted(:,1))); 
NoEpochs = length(FinalTTs);
nEst = 0;
nEst_c = 0;     % DGPS

for j = 1:NoEpochs
    for Iter = 1:MaxIter
        HTH = zeros(4,4);
        HTH_c = zeros(4,4);     % DGPS
        HTy = zeros(4,1);
        HTy_c = zeros(4,1);     % DGPS
        HT_ppcp = zeros(4,1);      % pp-CP
                
        indexQM = find(QM(:,1) == FinalTTs(j));
        QM_1 = QM(indexQM,:);
        indexPRC = find(PRC_Sorted(:,1) == FinalTTs(j));   % PRC 전체 행렬에서 현재 gs의 해 값 추출
        PRC_1 = PRC_Sorted(indexPRC,:);                    % PRC 전체 행렬에서 추출된 gs 행의 데이터 추출
        
        NoPRN(j,:) = [length(PRC_1(:,2)) length(QM_1(:,2))];
          %% prn 파일내에 obs prn이 없는 경우를 대비하기 위해서
        existprn = intersect(unique(PRC_1(:,2)), QM_1(:,2));
        arrSV = zeros(length(existprn),1);
        for k = 1:length(existprn)
            arrSV(k) = find(QM_1(:,2) == existprn(k,1));
        end
        QM_1 = QM_1(sort(arrSV),:);
        
        NoSats = length(QM_1);
        gs = QM_1(1,1);
              
        vec_site = x(1:3)';             
        vec_site_c = x_c(1:3)';     % DGPS
        ZHD = TropGPTh(vec_site, gw, gs);                 %: TROP: GPT
                
        for i = 1:NoSats
            prn = QM_1(i,2);
            obs = QM_1(i,4);   
            prc = PRC_1(find(PRC_1(:,2) == prn),3);   % DGPS
            icol = PickEPH(eph, prn, gs);
            toe = eph(icol, 8); a = eph(icol, 19); b = eph(icol, 20); c = eph(icol, 21); Tgd = eph(icol, 23);
            
            %----- 신호전달시간 계산
            STT = GetSTTbrdc(gs, prn, eph, vec_site);
            STT_c = GetSTTbrdc(gs, prn, eph, vec_site_c);    % DGPS
            tc = gs - STT;
            tc_c = gs - STT_c;                               % DGPS
                        
            %----- 위성궤도 계산
            vec_sat = GetSatPosNC(eph, icol, tc);          
            vec_sat_c = GetSatPosNC(eph, icol, tc_c);         % DGPS
            vec_sat = RotSatPos(vec_sat, STT);                %: 지구자전 고려  
            vec_sat_c = RotSatPos(vec_sat_c, STT_c);          %: DGPS 지구자전 고려
                                    
            %----- 최종 RHO 벡터 계산
            vec_rho = vec_sat - vec_site;
            vec_rho_c = vec_sat_c - vec_site_c;            % DGPS
            rho = norm(vec_rho);
            rho_c = norm(vec_rho_c);                          % DGPS
            [az,el] = xyz2azel(vec_rho, AppLat, AppLon);        
                        
            if el >= eleCut %15
                W = 1;
                dRel = GetRelBRDC(eph, icol, tc);
                dRel_c = GetRelBRDC(eph, icol, tc_c);                                           % DGPS
                
                dIono = ionoKlob(al, be, gs, az, el, vec_site);                                 % Klobuchar
                dTrop_G = ZHD2SHD(gw, gs, vec_site, el, ZHD);                                   % GPT
                
                dtSat = a + b*(tc - toe) + c*(tc - toe)^2 - Tgd + dRel;                         % relativistic
                dtSat_c = a + b*(tc_c - toe) + c*(tc_c - toe)^2 - Tgd + dRel_c;                 % DGPS
                
                com = rho + x(4) - CCC * dtSat + dIono + dTrop_G;                               % PP with correction
%                 com = rho + x(4) - CCC * dtSat                                                 % PP
                com_c = rho_c + x_c(4) - CCC * dtSat_c - prc;                                   % DGPS
                cp =  - dIono - dTrop_G - prc;                                                  % pp-CP_PRC with correction
%                 cp =  - prc;                                                                    % pp-CP_PRC
                
                y = obs - com;                                                                  % PP
                y_c = obs - com_c;                                                              % DGPS
                
                H = [ -vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 1];                       % PP
                H_c = [ -vec_rho_c(1)/rho_c -vec_rho_c(2)/rho_c -vec_rho_c(3)/rho_c 1];         % DGPS
                
                HTH = HTH + H'*W*H;                                                     % PP
                HTH_c = HTH_c + H_c'*W*H_c;                                             % DGPS
                
                HTy = HTy + H'*W*y;                                                     % PP
                HTy_c = HTy_c + H_c'*W*y_c;                                             % DGPS
                HT_ppcp = HT_ppcp + H'*W*cp;                                            % pp-CP_PRC
            end
                        
        end
        xhat = inv(HTH) * HTy;                                                          % PP
        xhat_c = inv(HTH_c) * HTy_c;                                                    % DGPS
        xhat_ppcp =  - inv(HTH) * HT_ppcp;                                % pp-CP
        x = x + xhat;                                                                   % PP
        x_c = x_c + xhat_c;                                                             % DGPS
        x_ppcp = x + xhat_ppcp;                                                         % pp-CP
        
%         if norm(xhat) < EpsStop;
        if Iter == 4;
            nEst = nEst + 1;
            estm(nEst,1) = gs;
            estm_c(nEst,1) = gs;
            estm_ppcp(nEst,1) = gs;
            estm(nEst,2:5) = x(1:4);
            estm_c(nEst,2:5) = x_c(1:4);
            estm_ppcp(nEst,2:5) = x_ppcp(1:4);
            break;
        end

    end
    
    user_gd(j,:) = xyz2gd(estm(j,2:4)); % user의 xyz값 gd 로 변환
    user_xyz(j,:) = estm(j,2:4);        % user xyz값 행렬로 변환 
  
end

%% 측위오차 분석 & 그래프 작성
estm = estm(1:nEst, :);
% [dXYZ, DGPSdXYZ, dNEV, DGPSdNEV] = PosTErrorsDGPS(estm(:, 1), TruePos, estm(:, 2:5), estm_c(:, 2:5));
[dXYZ, dNEV] = PosTErrors_result(estm(:, 1), TruePos, estm(:, 2:5));
[DGPSdXYZ, DGPSdNEV] = PosTErrors_result(estm_c(:, 1), TruePos, estm_c(:, 2:5));
[ppCPdXYZ, ppCPdNEV] = PosTErrors_result(estm_ppcp(:, 1), TruePos, estm_ppcp(:, 2:5));

% fid_out = fopen('estm15143.txt', 'w');
% fid_out_c = fopen('DGPSestm15143.txt', 'w');
% fid_out_ppcp = fopen('DGPSestm15143_ppCP.txt', 'w');
% for k = 1:nEst
%     fprintf(fid_out, '%8d %13.3f %13.3f %13.3f %8.3f \n', estm(k,: ));
%     fprintf(fid_out_c, '%8d %13.3f %13.3f %13.3f %8.3f \n', estm_c(k,: ));
%     fprintf(fid_out_ppcp, '%8d %13.3f %13.3f %13.3f %8.3f \n', estm_ppcp(k,: ));
% end
% fclose(fid_out);
% fclose(fid_out_c);
% fclose(fid_out_ppcp);

format long
[dNE dU d3] = RMS(dNEV); rmsdNEV = [dNE dU d3]
[DGPSdNE DGPSdU DGPSd3] = RMS(DGPSdNEV); rmsDGPSdNEV = [DGPSdNE DGPSdU DGPSd3]
[ppCPdNE ppCPdU ppCPd3] = RMS(ppCPdNEV); rmsppCPdNEV = [ppCPdNE ppCPdU ppCPd3]

% figure(98)
% hold on; grid on;
% axis square
% plot3(dNEV(:,2), dNEV(:,1), dNEV(:,3),'o'); 
% plot3(DGPSdNEV(:,2), DGPSdNEV(:,1), DGPSdNEV(:,3),'r.'); 
% plot3(ppCPdNEV(:,2), ppCPdNEV(:,1), ppCPdNEV(:,3),'m*'); 

figure(99)
hold on; grid on;
axis square
plot3(DGPSdNEV(:,2), DGPSdNEV(:,1), DGPSdNEV(:,3),'kd'...
    , ppCPdNEV(:,2), ppCPdNEV(:,1), ppCPdNEV(:,3),'m*'...
    ,dNEV(:,2), dNEV(:,1), dNEV(:,3),'o'); 
legend('DGPS', 'PP-CP','PointPosition')
xlabel({['  rms DGPS \Delta NE = ', num2str(rmsDGPSdNEV(1))...
    ,'   rms DGPS \Delta V = ', num2str(rmsDGPSdNEV(2))...
    ,'   rms DGPS \Delta 3D = ', num2str(rmsDGPSdNEV(3))],...
    ['   rms PP-CP \Delta NE = ', num2str(rmsppCPdNEV(1))...
    ,'   rms PP-CP \Delta V = ', num2str(rmsppCPdNEV(2))...
    ,'   rms PP-CP \Delta 3D = ', num2str(rmsppCPdNEV(3))],...
        ['   rms PP \Delta NE = ', num2str(rmsdNEV(1))...
    ,'   rms PP \Delta V = ', num2str(rmsdNEV(2))...
    ,'   rms PP \Delta 3D = ', num2str(rmsdNEV(3))]});