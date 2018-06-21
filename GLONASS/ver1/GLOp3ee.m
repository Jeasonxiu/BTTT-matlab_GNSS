%% GLONASS 방송궤도력 PPP 코드정리 %% 2015.01.16
clear all; tic;
%% 상수 정의
CCC = 299792458.;   % CCC = Speed of Light [m/s]
%% QM 파일, 날짜 핸들링
DOY = 212; YY = 14;
% DOY = 225; YY = 16;
% --- GPS WEEK 결정
[gw, gd] = ydoy2gwgd(YY, DOY);
% --- QM file 준비
FileQM = 'IHUR2120';
% FileQM = 'uatA';
[arrQM, FinalPRNs, FinalTTs] = ReadQM(FileQM);

ObsType = 220; % 관측치 설정: GLONASS[211(L1), 212(L2), 220(C1), 221(P1), 222(P2), 231(D1), 232(D2)]
% ObsType = 320; % 관측치 설정: GLONASS[211(L1), 212(L2), 220(C1), 221(P1), 222(P2), 231(D1), 232(D2)]

QM = SelectQM(arrQM,ObsType);
FinalPRNs = unique(QM(:,2));
FinalTTs = unique(QM(:,1));
% --- 좌표 입력
TruePos = [-3026676.0349  4067187.8095  3857246.8615 ];
% TruePos = [-3041235.578 4053941.677 3859881.013];   % : JPspace A point
%% 사용 파일 입력 및 데이터 추출
% --- GLONASS 방송궤도력 파일: EPH, LeapSeceond, Tau_C
FileNavGLO = strcat('brdc', num2str(DOY), '0.', num2str(YY), 'g');
EphGlo = ReadEPH_glo(FileNavGLO);
Tau_c = ReadTauC(FileNavGLO);
% LeapSec = GetLeapSec(FileNavGLO);
% --- GPS 방송궤도력 파일: 전리층-Klobuchar 변수
FileNav = strcat('brdc', num2str(DOY), '0.', num2str(YY), 'n');
[al, be] = GetALBE(FileNav);
LeapSec = GetLeapSec(FileNav);
% --- IONEX 파일: DCB
FileIon = strcat('igsg', num2str(DOY), '0.', num2str(YY), 'i');
DCB = ReadDCB(FileIon);

%% 좌표 초기치 설정 및 위경도 변환
AppPos = TruePos;
gd = xyz2gd(AppPos); AppLat = gd(1); AppLon = gd(2);
%% 추정을 위한 매개변수 설정
Maxiter = 5;
EpsStop = 1e-4;
ctr = 1; eleCut=15;
x = [AppPos ctr]; x = x'; deltat = 15;
%%
NoEpochs = length(FinalTTs);
estm = zeros(NoEpochs, 5);
nEst = 0;

for j = 1:NoEpochs
    for iter = 1:Maxiter
        HTH = zeros(4,4);
        HTy = zeros(4,1);
        
        gs = FinalTTs(j);
        indexQM = find(QM(:,1) == gs);
        QM1e = QM(indexQM,:);
        NoSats = length(QM1e);
        
        vec_site = x(1:3)';
        % --- 연직방향 대류권 지연량 계산; 에포크별 한 번 계산
        ZHD = TropGPTh(vec_site, gw, gs); % GPT
        
        for i = 1:NoSats
            prn = QM1e(i,2);
            obs = QM1e(i,4);
            icol = PickEPH_GLO(EphGlo, prn, gs);
            
            TauN = EphGlo(icol,12); % : tau & gamma - 시계오차 보정에 필요
            GammaN = EphGlo(icol,13);
            ch_num = EphGlo(icol,16); % : channel number - 전리층 보정에 필요
            
            %%
            STT = GetSTTbrdc_GLO(gs, EphGlo, icol, x(1:3), deltat); % 신호전달시간 계산
            tc = gs - STT; % 신호전달시간 보정
            tc = tc - LeapSec + Tau_c; % : LeapSecond & TauC 보정
            
            % 신호전달시간 & LeapSecond 보정한 위성 위치 산출
            [sat_pos,sat_vel] = GetSatPosGLO(EphGlo,icol,tc,deltat);
            SatPos = sat_pos;
%             SatPos = PZ2WGS(sat_pos);
            SatPos = RotSatPos(SatPos,STT); % 지구 자전효과 고려
%             SatPos = SatPos';
            
            DistXYZ = SatPos - vec_site;
            DistNorm = norm(DistXYZ);
            
            [az,el] = xyz2azel(DistXYZ, AppLat, AppLon);
            %             if el>=eleCut
            %%
            dIono = Klo_R(vec_site, al, be, tc, SatPos,ch_num); % 전리층 보정
            dTrop = ZHD2SHD(gw, gs, vec_site, el, ZHD); % 대류권 보정
            
            SatVel = sat_vel;
            dRel = (-2/CCC^2) * dot(SatPos, SatVel); % 상대론적 효과
            dDCB = AppDCB_glo(DCB,prn);% DCB
            %%
            tsv = tc;
            tb = EphGlo(icol,2) + LeapSec; % GPStime - 16; / GLOtime + 16; 확인하기
            dtSat = TauN - GammaN*(tsv-tb) + Tau_c + dRel + dDCB;
            com = DistNorm + x(4) - CCC*dtSat + dTrop + dIono;
            %%
            y = obs - com;
            
            H = [-DistXYZ(1)/DistNorm -DistXYZ(2)/DistNorm -DistXYZ(3)/DistNorm 1];
            HTH = HTH + H'*H;
            HTy = HTy + H'*y;
            %             end  % if el>=eleCut
        end % --- for i = 1:NoSats
        xhat = inv(HTH) * HTy;
        x = x + xhat;
        if norm(xhat) < EpsStop;
            nEst = nEst + 1;
            estm(nEst,1) = gs;
            estm(nEst,2:5) = x(1:4);
            fprintf('%8d : %8.2f\n', j, norm(x(1:3)' - TruePos))
            break;
        end
    end % --- for iter = 1:Maxiter
end % % --- for j = 1:NoEpochs

% figure(200)
% plot(estm(:,1),estm(:,5),'-ob'); grid on;
toc; 

%% Analysis of Positioning Errors
estm = estm(1:nEst, :);
% [dXYZ, dNEV]=PosErrors(estm(:,1), TruePos,estm(:,2:4));
[dXYZ, dNEV] = PosTErrors2(estm(:,1), TruePos, estm(:,2:5));