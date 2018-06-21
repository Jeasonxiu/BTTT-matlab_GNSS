function [] = PP(obsfile, navfile)
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
%% 불변 변수 설정: 빛의 속도, 관측치
CCC = 299792458.;   % CCC = Speed of Light [m/s]
ObsType = 120;      % 사용할 관측치 설정 - 120: C/A = C1
%% 임계고도각 설정
eleCut = 15;

%% rinex file 입력
% obsfile = 'ubx1_150707_10m_ubx.obs';
% navfile = 'ubx1_150707_10m_ubx.nav';
%% Observation Rinex 파일로 부터 QMfile 생성
WriteObs(obsfile)

%% 만들어진 QMfile을 obsevation Rinex 파일 관측 날짜, 시간으로 변경
rename = renameQMfile(obsfile,'HYU');

%% QM 파일 읽어들여서 행렬로 저장하고, 사용할 관측치 추출
[arrQM, FinalPRNs, FinalTTs] = ReadQM(rename); 
QM = SelectQM(arrQM, ObsType);

%% 항법메시지를 읽어들여서 행렬로 저장하고, Klobuchar 모델 추출
eph = ReadEPH(navfile);

%% 라이넥스 파일에서 대략적인 관측소 좌표를 뽑아내고 위경도로 변환
AppPos = GetAppPos(obsfile);
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
        NoSats = length(QM_1);
        gs = QM_1(1,1);
              
        vec_site = x(1:3)';
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
            vec_rho = (vec_sat - vec_site)';
            rho = norm(vec_rho);
            [az,el] = xyz2azel(vec_rho, AppLat, AppLon);        
            if el >= eleCut %15
                W = 1;
                dRel = GetRelBRDC(eph, icol, tc);
                dtSat = a + b*(tc - toe) + c*(tc - toe)^2 - Tgd + dRel;
                com = rho + x(4) - CCC * dtSat;
                y = obs - com;
                H = [ -vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 1];
                HTH = HTH + H'*W*H;
                HTy = HTy + H'*W*y;
                
            end
        end
        xhat = inv(HTH) * HTy;
        x = x + xhat;
        
        if norm(xhat) < EpsStop;
            nEst = nEst + 1;
            estm(nEst,1) = gs;
            estm(nEst,2:5) = x(1:4);
%             fprintf('%8d : %8.2f\n',j, norm(x(1:3)' - TruePos))
            break;
        end
        
    end
    
    user_gd(j,:) = xyz2gd(estm(j,2:4)); % user의 xyz값 gd 로 변환
    user_xyz(j,:) = estm(j,2:4);        % user xyz값 행렬로 변환 
end

%% user_position mean값 기준 topology 계산 과정
[user_mean] = PlotMeanTopo(user_xyz);

figure(201)
grid on
hold on
plot(user_gd(:,2), user_gd(:,1),'ro','markeredgecolor','y','markerfacecolor','r','markersize',3)
plot(mean(user_gd(:,2)), mean(user_gd(:,1)), 'b*','markeredgecolor','b','markerfacecolor','r','markersize',5)
plot_google_map;


