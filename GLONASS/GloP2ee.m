%%% GLONASS 방송궤도력 PPP 코드정리 %%%
% August 24th, 2016, Joonseong Gim, 기존측위코드에서 RK4, PZ90.02 -> PZ90.11 수정


%% GLONASS 방송궤도력 PPP 코드정리 %% 2015.01.16
clear all; tic;
%% 상수 정의
CCC = 299792458.;   % CCC = Speed of Light [m/s]
%% QM 파일, 날짜 핸들링
% DOY = 212; YY = 14;
% DOY = 225; YY = 16;
% DOY = 237; YY = 16;
DOY = 153; YY = 16;
% --- GPS WEEK 결정
[gw, GD] = ydoy2gwgd(YY, DOY);
% --- QM file 준비
FileQM = 'IHUR2120';
% FileQM = 'uatA';
% FileQM = 'QSDT1_16_ob082';
FileQM = 'QDRUB_ob153_3';
[arrQM, FinalPRNs, FinalTTs] = ReadQM(FileQM);
ObsType = 320; % 관측치 설정: GLONASS[211(L1), 212(L2), 220(C1), 221(P1), 222(P2), 231(D1), 232(D2)]
QM = SelectQM(arrQM,ObsType);
QM2 = SelectQM(arrQM,341);
FinalPRNs = unique(QM(:,2));
FinalTTs = unique(QM(:,1));
% --- 좌표 입력
% TruePos = [-3026676.0349  4067187.8095  3857246.8615 ];
TruePos = [-3041235.578 4053941.677 3859881.013];   % : JPspace A point
% TruePos = [-3041241.741 4053944.143 3859873.640];   % : JPspace B point
AppPos = TruePos;
% --- 좌표 위경도 변환
gd = xyz2gd(AppPos); AppLat = gd(1); AppLon = gd(2);
%% 사용 파일 입력 및 데이터 추출
% --- GLONASS 방송궤도력 파일 : EphGLO, LeapSecond, TauC
FileNavGLO = strcat('brdc',num2str(DOY),'0.',num2str(YY),'g');
EphGlo = ReadEPH_glo(FileNavGLO);
TauC = ReadTauC(FileNavGLO); %%
% --- GPS 방송궤도력 파일: 전리층-Klobuchar
FileNav = strcat('brdc',num2str(DOY),'0.',num2str(YY),'n');
[al,be] = GetALBE(FileNav); %%
LeapSec = GetLeapSec(FileNav); %%
% --- IONEX 파일: DCB
FileIon = strcat('igsg',num2str(DOY),'0.',num2str(YY),'i');
DCB = ReadDCB(FileIon);
%% 추정을 위한 매개변수 설정
Maxiter = 10; 
EpsStop = 1e-5;
ctr = 1; eleCut=15; deltat = 1;
x = [AppPos ctr]; x = x';

%%
NoEpochs = length(FinalTTs);
estm = zeros(NoEpochs,5);
nEst = 0;
for j = 1:NoEpochs
% for j = 1:85
    indexQM = find(QM(:,1) == FinalTTs(j));
    QM_1 = QM(indexQM,:); NoSats = length(QM_1);
    QM_2 = QM2(indexQM,:); NoSats = length(QM_2);
    gs = QM_1(1,1);
    vec_site = x(1:3)';
    visiSat(j,1) = gs; visiSat(j,2) = NoSats;
    GloQM_snr = QM_2(find(QM_2(:,3) == 341),:); GloSats = length(GloQM_snr(:,4));     % GLO C1
    %% 1초 간격으로 정각, 30분에 위성 위치 계산(방송궤도력 15분, 45분/Forward,Backward)
    sT= mod(gs,86400)/3600;
    sTh = floor(sT); sTm = sT - sTh;
    if j == 1 % 시작시간 위성 위치 array
        [Sat_ar] = GetSatPosGLO_new(EphGlo,gs,deltat);
%         [Sat_ar] = GetSatPosGLO_my(EphGlo,gs,deltat);
        %         fprintf('시작 시간 위치 계산\n');
    elseif (j ~= 1 && sTm == 0) % 정각일 때
        clear Sat_ar
        [Sat_ar] = GetSatPosGLO_new(EphGlo,gs,deltat);
%         [Sat_ar] = GetSatPosGLO_my(EphGlo,gs,deltat);
        %         fprintf('%d시 정각 갱신\n',sTh);
    elseif (j ~= 1 && mod(sTm,0.5) == 0) % 30분 일 때
        clear Sat_ar
        [Sat_ar] = GetSatPosGLO_new(EphGlo,gs,deltat);
%         [Sat_ar] = GetSatPosGLO_my(EphGlo,gs,deltat);
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
            S1 = GloQM_snr(find(GloQM_snr(:,2) == prn),4);
            tc = gs - LeapSec;
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
            if el>=eleCut
                % 전리층 보정 % 대류권 보정
                %                 W = MakeW_elsnr(el,S1);
                %                 W = MakeW_elpr(el);
                W = 1;
                ttc = tc - TauC;
                dIono = Klo_R(vec_site,al,be,ttc,SatPos,ch_num);
                dTrop = ZHD2SHD(gw,gs,vec_site,el,ZHD);
                % 상대적 효과 (→위성속도 이용)
                dRel = (-2/CCC^2) * dot(SatPos, SatVel);
                % DCB 고려
                dDCB = AppDCB_glo(DCB,prn+100);
                %%
                % 위성시계오차 tsv, tb
                tsv = tc; tb = EphGlo(icol,2) + LeapSec; % GPStime - 16; / GLOtime + 16; 확인하기
%                 dtSat = TauN - GammaN*(tsv-tb) + TauC + dRel + dDCB;
                dtSat = TauN - GammaN*(tsv-tb) + TauC + dRel;
                
                com = DistNorm + x(4) - CCC*dtSat + dTrop + dIono;
                y = obs - com;
                H = [-DistXYZ(1)/DistNorm -DistXYZ(2)/DistNorm -DistXYZ(3)/DistNorm 1];
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
            fprintf('gs: %6.0f     %2.5f \n',gs,norm(TruePos'-x(1:3)));
            break;
        end
    end
    user_gd(j,:) = xyz2gd(estm(j,2:4)); % user의 xyz값 gd 로 변환
end
estm = estm(find(estm(:,1) > 0),:);
for es = 1:length(estm(:,1))
    user_gd(es,:) = xyz2gd(estm(es,2:4)); % user의 xyz값 gd 로 변환
    AppLat = user_gd(es,1); AppLon = user_gd(es,2);
    user_xyz(es,:) = estm(es,2:4);        % user xyz값 행렬로 변환
end
%% Analysis of Positioning Errors
% estm = estm(1:nEst, :);
% [dXYZ, dNEV]=PosErrors1(EstPos(:,1), TruePos, EstPos(:,2:4),visiSat(:,1));


%% 구글 plot
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




%% PostErrors
% [dXYZ, dNEV] = PosTErrors5(estm, UBLOX, visiSat);
% 
% vrsfile = '160524_1_adm.txt';
% VRS = load(vrsfile);
% VRS(:,1) = VRS(:,1) +17;
% UBLOX(:,1) = VRS(:,1);
% UBLOX(:,2:4) = VRS(:,5:7);
% [dXYZ, dNEV] = PosTErrors5(estm, UBLOX, visiSat);
% [dXYZ, dNEV] = PosTErrors5(estm_d, UBLOX, visiSat);
% 
% VRS_text = 'SDT1_VRS_16237.txt';
% [target,UBLOX] = gapconv(VRS_text, 0.43, 16, 237);


%% PostErrors
[dXYZ, dNEV]=PosTErrors4(estm(:,1), TruePos, estm(:,2:5),visiSat);

%% 위성별 SNR
% PPPlotQM(QM2,341)

toc;