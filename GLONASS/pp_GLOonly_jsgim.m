%% Clearing everything before each run
clear all; close all;
%% 기초 설정: 빛의 속도, 관측치
CCC = 299792458.;   %: CCC = Speed of Light [m/s]
ObsType = 320;      %: GLONASS C1
LeapSec = 18;
eleCut  = 5;
%% 날짜와 사이트 
YY = 17; DOY = 025; QMfile = 'QM170125_A'; TruePos = [-3041235.578 4053941.677 3859881.013];       % 대성 A point
%% 
FileNav = strcat('brdm', num2str(DOY,'%03d'), '0.', num2str(YY,'%02d'), 'p');   %: Navigation RINEX file
%% QM과 eph_glo는 저장해두자
% [arrQM, dummy, dummy] = ReadQM(QMfile);
% [QM, FinalPRNs, FinalTTs] = SelectQMFinal(arrQM, ObsType);
% [eph, trashPRN, trashT] = ReadEPH_all(FileNav);
% ephGLO = eph;
% ephGLO = ephGLO(ephGLO(:,1) < 400 & ephGLO(:,1) > 300,:);
load 'ws_ppGLO';
%% 
[gw, gd] = ydoy2gwgd(YY, DOY); %: GPS WEEK 결정
TauC = ReadTauC2(FileNav);
%% 좌표 초기치 설정 및 위경도 변환
AppPos = TruePos;
gds = xyz2gd(AppPos); AppLat = gds(1); AppLon = gds(2);
%% 추정을 위한 매개변수 설정
Maxiter = 5;
EpsStop = 1e-4;
cdtr = 1; 
deltat = 60;
x = [AppPos cdtr]';         

NoEpochs = length(FinalTTs);
estm = zeros(NoEpochs, 5);
%%  이전 epoch과 현재 epoch을 비교해서, icol의 변화가 있는지 판단하기 위한 arr생성 파트 입니다.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SatPosArr_before = zeros(24,25); %:?? 이줄과 다음줄이 뭐하는 역할?
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
            
            icol = PickEPH_GLO2(ephGLO, prn+300, gs); %: PRN번호 세팅
            icolArr(i,1) = prn+300; 
            icolArr(i,2) = icol;
            TauN = ephGLO(icol,12); %: tau & gamma - 시계오차 보정에 필요
            GammaN = ephGLO(icol,13);
            jcol = find(SatPosArr_before(:,1) == prn+300); %:PRN번호 세팅
            if isempty(jcol); jcol = 0; end
            jcol = jcol(1);
            %% 위성위치 계산 파트
            if icol == icolArr_before(icolArr_before(:,1)==prn+300,2) % icol 변화 없으면 이전 epoch으로 계산
                STT = 0.075;   
                tc = gs - STT;
                tc = tc - LeapSec - TauC;
                [SatPos,SatVel,SatLS] = GetSatPosGLO_ver3(SatPosArr_before,jcol,tc,deltat); % 이전 epoch값으로 위성위치 계산
                SatPos = RotSatPos(SatPos, STT); % 지구자전효과 고려
            else % icol 변화시(새위성 또는 새로운 방송궤도력) 방송궤도력으로 계산
                STT = 0.075;
                tc = gs - STT;
                tc = tc - LeapSec - TauC;
                [SatPos,SatVel,SatLS] = GetSatPosGLO(ephGLO,icol,tc,deltat); % 방송궤도력으로 위성위치 계산
                SatPos = RotSatPos(SatPos, STT); % 지구자전효과 고려
            end
            icolArr_before(i,:) = icolArr(i,:);
            SatPosArr_before(i,:) = ephGLO(icol,:);
            
            vec_rho = SatPos - vec_site';
            rho = norm(vec_rho);
            [az,el] = xyz2azel(vec_rho, AppLat, AppLon);
            dIono = 0;
            dTrop = ZHD2SHD(gw, gs, TruePos, el, ZHD); % 대류권 보정
            dRel = (-2/CCC^2) * dot(SatPos, SatVel); % 상대론적 효과
            tb = ephGLO(icol,2) + LeapSec;
            dtSat = TauN - GammaN*(tc-tb) + TauC + dRel;
            
            if el >= eleCut % 임계고도각
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
                estm(nEst,6:7) = [NoSats, NoSatsUsed];      % 관측 위성 수와 사용 위성
                estm(nEst,8:10) = xyz2gd(x(1:3));           % 계산된 ecef를 gd로 변환
                fprintf('%8d : %3d : %8.2f : %8.2f : %8.2f\n', j, iter, x(1)' - TruePos(1), x(2)' - TruePos(2), x(3)' - TruePos(3));
                break;
            end
        else
            break;
        end
    end
end
%% 측위오차 분석 & 그래프 작성
estm = estm(1:nEst, :);
[dXYZ, dNEV] = PosTErrors2(estm(:,1), TruePos, estm(:,2:5),estm(:,6:7));

% %% 구글지도plot
% figure(99)
% hold on; grid on;
% plot(estm(:,10),estm(:,9),'bo')
% plot_google_map;