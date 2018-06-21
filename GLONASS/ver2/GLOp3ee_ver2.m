%%% GLONASS 방송궤도력 PPP 코드정리 %%%
% April 10th, 2015, Mi-So Kim, 기존측위코드에서 속도 개선함

clc; clear all; tic;
%% 상수 정의
CCC = 299792458.;   % CCC = Speed of Light [m/s]
%% QM 파일, 날짜 핸들링
DOY = 212; YY = 14;
% --- GPS WEEK 결정
[gw, gd] = ydoy2gwgd(YY,DOY); %%
% --- QM File 준비
FileQM = 'IHUR2120';
[arrQM, FinalPRNs, FinalTTs] = ReadQM(FileQM );
ObsType = 220; % 관측치 설정: [GLONASS] 211(L1), 212(L2), 220(C1), 221(P1), 222(P2),231(D1), 232(D2)
QM = SelectQM(arrQM,ObsType);
FinalPRNs = unique(QM(:,2));
FinalTTs = unique(QM(:,1));
%% 좌표 초기치 설정 및 위경도 변환
% --- 좌표 입력
TruePos = [-3026676.0349  4067187.8095  3857246.8615]; %IHUR1980.14o
AppPos = TruePos;
% --- 좌표 위경도 변환
gd = xyz2gd(AppPos); AppLat = gd(1); AppLon = gd(2);
%% 사용 파일 입력 및 데이터 추출
% --- GLONASS 방송궤도력 파일 : EphGLO, LeapSecond, TauC
FileNavGLO = strcat('brdc',num2str(DOY),'0.',num2str(YY),'g');
EphGlo = ReadEPH_GLO(FileNavGLO);
TauC = ReadTauC(FileNavGLO); %%
% --- GPS 방송궤도력 파일: 전리층-Klobuchar
FileNav = strcat('brdc',num2str(DOY),'0.',num2str(YY),'n');
[al,be] = GetALBE(FileNav); %%
LeapSec = GetLeapSec(FileNav); %%
% --- IONEX 파일: DCB
FileIon = strcat('igsg',num2str(DOY),'0.',num2str(YY),'i');
DCB = ReadDCB(FileIon);
%% 추정을 위한 매개변수 설정
Maxiter = 10; EpsStop = 1e-5;
ctr = 1; eleCut=15; deltat = 1;
x = [AppPos ctr]; x = x';

%%
NoEpochs = length(FinalTTs);
EstPos = zeros(NoEpochs,5);
nEst = 0;
for j = 1:NoEpochs
    indexQM = find(QM(:,1) == FinalTTs(j));
    QM_1 = QM(indexQM,:); NoSats = length(QM_1);
    gs = QM_1(1,1); 
    vec_site = x(1:3)';
    visiSat(j,2) = NoSats; visiSat(j,1) = gs;
    %% 1초 간격으로 정각, 30분에 위성 위치 계산(방송궤도력 15분, 45분/Forward,Backward)
    sT= mod(gs,86400)/3600;
    sTh = floor(sT); sTm = sT - sTh;    
    if j == 1 % 시작시간 위성 위치 array
%         [Sat_ar] = GetSatPosGLO_new(EphGlo,gs,deltat);
        [Sat_ar] = GetSatPosGLO_my(EphGlo,gs,deltat);
%         fprintf('시작 시간 위치 계산\n');
    elseif (j ~= 1 && sTm == 0) % 정각일 때
        clear Sat_ar
%         [Sat_ar] = GetSatPosGLO_new(EphGlo,gs,deltat);
        [Sat_ar] = GetSatPosGLO_my(EphGlo,gs,deltat);
%         fprintf('%d시 정각 갱신\n',sTh);
    elseif (j ~= 1 && mod(sTm,0.5) == 0) % 30분 일 때
        clear Sat_ar
%         [Sat_ar] = GetSatPosGLO_new(EphGlo,gs,deltat);
        [Sat_ar] = GetSatPosGLO_my(EphGlo,gs,deltat);
%         fprintf('%d시 30분 갱신\n',sTh);
    end
    % --- 연직방향 대류권 지연략 계산: 에포크별 한 번 계산
    ZHD = TropGPTh(vec_site, gw, gs); % GPT
    %% 
    for iter = 1:Maxiter
        HTH = zeros(4,4);
        HTy = zeros(4,1);
        for i = 1:NoSats
            prn = QM_1(i,2);
            obs = QM_1(i,4);
            
            tc = gs - 16;
            icol=PickEPH_GLO2(EphGlo, prn, tc);
            
            TauN=EphGlo(icol,12); GammaN=EphGlo(icol,13); %: tau & gamma 시계오차 보정에 사용
            ch_num=EphGlo(icol,16); %: channel number 전리층 보정에 사용
            
            % 신호전달시간 계산
            STT = GetSTTbrdcGLO2(Sat_ar,gs,prn,x(1:3));           
            % LeapSecond & 신호전달 시간을 보정한 위성 위치 산출
            [SatPos, SatVel] = SatPosLS_STT(Sat_ar,gs,prn,LeapSec,STT,TauC);

            % 지구 자전효과 고려
            SatPos = RotSatPos(SatPos,STT);
%             SatPos = SatPos';
            
            DistXYZ = SatPos - x(1:3)';
            DistNorm = norm(DistXYZ);           
            [az,el] = xyz2azel(DistXYZ, AppLat, AppLon);
            %%
%             if el>=eleCut
                % 전리층 보정 % 대류권 보정
                ttc = tc - TauC;
                dIono = Klo_R(vec_site,al,be,ttc,SatPos,ch_num); 
                dTrop = ZHD2SHD(gw,gs,vec_site,el,ZHD);
                % 상대적 효과 (→위성속도 이용)
                dRel = (-2/CCC^2) * dot(SatPos, SatVel);
                % DCB 고려
                dDCB = AppDCB_glo(DCB,prn);
                %%
                % 위성시계오차 tsv, tb
                tsv = tc; tb = EphGlo(icol,2) + 16; % GPStime - 16; / GLOtime + 16; 확인하기
%                 dtSat = TauN - GammaN*(tsv-tb) + TauC + dRel + dDCB;
                dtSat = TauN - GammaN*(tsv-tb) + TauC + dRel ;
                
                com = DistNorm + x(4) - CCC*dtSat + dTrop + dIono;
                y = obs - com;
                H = [-DistXYZ(1)/DistNorm -DistXYZ(2)/DistNorm -DistXYZ(3)/DistNorm 1];
                HTH = HTH + H'*H;
                HTy = HTy + H'*y;
%             end
        end
        xhat = inv(HTH) * HTy;
        x = x + xhat;
        if norm(xhat) < EpsStop;
            nEst = nEst + 1;
            EstPos(nEst,1) = gs;
            EstPos(nEst,2:5) = x(1:4);
            fprintf('gs: %6.0f     %2.5f \n',gs,norm(TruePos'-x(1:3)));    
            break;
        end
    end
end

%% Analysis of Positioning Errors
EstPos = EstPos(1:nEst, :);
[dXYZ, dNEV]=PosErrors1(EstPos(:,1), TruePos, EstPos(:,2:4),visiSat(:,2));
toc;
[dXYZ, dNEV]=PosTErrors3(EstPos(:,1), TruePos, EstPos(:,2:5),visiSat);
