% function [estm] = PP(obsfile, navfile)
%
%function [] = PP(obsfile, navfile)
%
%   Read the files(obs, nav) and Plot 'horizonal error', 'Vertical
%   error', from user_mean position(x,y,z)
%   
%   input obsfile : Rinex Observation file
%   input navfile : Rinex navigation file
%   
%   Example : PP('*.15n', '*.15n')
%
%   coded by Joonseong Gim, Jan 13, 2016
%

clear all;


%% 불변 변수 설정: 빛의 속도, 관측치
CCC = 299792458.;   % CCC = Speed of Light [m/s]
ObsType = 120;      % 사용할 관측치 설정 - 120: C/A = C1

% TruePos = [-3041235.578 4053941.677 3859881.013];   % : JPspace A point
TruePos = [-3041241.741 4053944.143 3859873.640];   % : JPspace B point
%% 임계고도각 설정
eleCut = 15;

%% rinex file 입력
obsfile = 'DBUU049i_16.16o';
navfile = 'brdc0490.16n';
% obsfile = 'DAEB011c.16o';
% navfile = 'brdc0110.16n';
%% Observation Rinex 파일로 부터 QMfile 생성
WriteObs(obsfile)

%% 만들어진 QMfile을 obsevation Rinex 파일 관측 날짜, 시간으로 변경
rename = renameQMfile(obsfile);
[YY, DOY] = obs2YYDOY(obsfile);
[gw, GD] = ydoy2gwgd(YY, DOY); %: GPS WEEK 결정
%% QM 파일 읽어들여서 행렬로 저장하고, 사용할 관측치 추출
[arrQM, FinalPRNs, FinalTTs] = ReadQM(rename); 
QM = SelectQM(arrQM, ObsType);
%% PRC load
PRCfile = 'JPRT160218.t1';
[PRC_Sorted] = PRCsort(PRCfile, QM); 

%% 항법메시지를 읽어들여서 행렬로 저장하고, Klobuchar 모델 추출
eph = ReadEPH(navfile);
[al, be] = GetALBE(navfile);
%% 라이넥스 파일에서 대략적인 관측소 좌표를 뽑아내고 위경도로 변환
AppPos = GetAppPos(obsfile);
if AppPos(1) == 0
    AppPos = TruePos;
end

gd = xyz2gd(AppPos); AppLat = gd(1); AppLon = gd(2); 

%% 추정에 필요한 초기치 설정
MaxIter = 4;
EpsStop = 1e-5;
ctr = 1;
x = [AppPos ctr]; x = x';

%% 추정과정 시작
NoEpochs = length(FinalTTs);
nEst = 0;


for j = 1:NoEpochs
    for Iter = 1:MaxIter
        HTH = zeros(4,4);
        HTy = zeros(4,1);
        
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
        ZHD = TropGPTh(vec_site, gw, gs);                 %: TROP: GPT
        
        for i = 1:NoSats
            prn = QM_1(i,2);
            obs = QM_1(i,4);      
            icol = PickEPH(eph, prn, gs);
            toe = eph(icol, 8); a = eph(icol, 19); b = eph(icol, 20); c = eph(icol, 21); Tgd = eph(icol, 23);
            %----- 신호전달시간 계산
            STT = GetSTTbrdc(gs, prn, eph, vec_site);
            tc = gs - STT;
            %----- 위성궤도 계산
            vec_sat = GetSatPosNC(eph, icol, tc);          
            vec_sat = RotSatPos(vec_sat, STT);                      %: 지구자전 고려  
            %----- 최종 RHO 벡터 계산
            vec_rho = vec_sat - vec_site;
            rho = norm(vec_rho);
            [az,el] = xyz2azel(vec_rho, AppLat, AppLon); 
            
            if el >= eleCut %15
                W = 1;
                dRel = GetRelBRDC(eph, icol, tc);
                dtSat = a + b*(tc - toe) + c*(tc - toe)^2 - Tgd + dRel;
                dIono = ionoKlob(al, be, gs, az, el, vec_site);
%                 dTrop_H = Hopfield(el, 11, vec_site, 9999);                   % Hopfield Model
%                 dTrop_S = Saastamoinen(el, 11, vec_site, 9999);               % Saastamoinen Model
                dTrop_G = ZHD2SHD(gw, gs, vec_site, el, ZHD);                   % GPT model
%                 com = rho + x(4) - CCC * dtSat + dIono + dTrop_H;             % Hopfield Model
%                 com = rho + x(4) - CCC * dtSat + dIono + dTrop_S;             % Saastamoinen Model
                com = rho + x(4) - CCC * dtSat + dIono + dTrop_G;             % GPT Model
%                 com = rho + x(4) - CCC * dtSat;
                y = obs - com;
                H = [ -vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 1];
                HTH = HTH + H'*W*H;
                HTy = HTy + H'*W*y;
                
            end
        end
        xhat = inv(HTH) * HTy;
        x = x + xhat;
        
        if Iter == 4;
            nEst = nEst + 1;
            estm(nEst,1) = gs;
            estm(nEst,2:5) = x(1:4);
%             fprintf('%8d : %8.2f\n',j, norm(x(1:3)' - TruePos))
            break;
        end
        
    end
    
    user_gd(j,:) = xyz2gd(estm(j,2:4)); % user의 xyz값 gd 로 변환
    AppLat = user_gd(j,1); AppLon = user_gd(j,2);
    user_xyz(j,:) = estm(j,2:4);        % user xyz값 행렬로 변환 
end

%% user_position mean값 기준 topology 계산 과정
% [user_mean] = PlotMeanTopo(user_xyz);
% user_mean = [mean(user_xyz(:,1)) mean(user_xyz(:,2)) mean(user_xyz(:,3))];
user_mean = [-3050139.16635833,4040764.76856110,3866462.84760506];

% [dXYZ, dNEV] = PosTErrorsJOON(estm(:, 1), user_mean, estm(:, 2:5));
[dXYZ, dNEV] = PosTErrorsJOON(estm(:, 1), TruePos, estm(:, 2:5));


% figure(201)
% axis([127.03 127.041 37.528 37.548]);
% plot_google_map;
% axis equal
% axis([127.03 127.041 37.528 37.548]);
% 
% for i = 1:length(user_gd)
%     GD_lon(i,1) = user_gd(i,2);
%     GD_la(i,1) = user_gd(i,1);
%     figure(201)
%     grid on
%     hold on
%     plot(GD_lon, GD_la,'ro','markeredgecolor','y','markerfacecolor','r','markersize',3)
%     drawnow
% end

% figure(201)
% grid on
% hold on
% % axis([(min(user_gd(:,2))-0.0001) (max(user_gd(:,2))+0.0001) (min(user_gd(:,1))-0.0001) (max(user_gd(:,1))+0.0001)]);
% plot(user_gd(:,2), user_gd(:,1),'ro','markeredgecolor','y','markerfacecolor','r','markersize',3)
% plot_google_map;