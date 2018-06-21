%% �ڵ��ǻ�Ÿ� �������� �˰���
% 07/01/2016 : Joonseong
close all; clear all;

%% �Һ� ���� ���� : ���� �ӵ�, ����ġ
CCC = 299792458.;   % CCC = Speed of Light [m/s]
ObsType = 120;      % ����� ������ ���� - 120 : C/A = C1
ObsType2 = 141;      % ����� ����ġ ���� - 141: snr = S1

%% �Ӱ���� ����
eleCut = 15;

%% QM ���� �ڵ鸵
% obsfileBs = 'DDBs115_31.16o';                                 % Reference Obsevation
% navfileBs = 'DDBs115_31.16n';                                 % Reference Navigation
% obsfileRv = 'DDRv115_11.16o';                                 % estm Obsevation
% navfileRv = 'DDRv115_11.16n';                                 % estm Navigation
obsfileBs = 'SBBs055_14.16o';                                 % Reference Obsevation
navfileBs = 'SBBs055_14.16n';                                 % Reference Navigation
obsfileRv = 'SBRv055_34.16o';                                 % estm Obsevation
navfileRv = 'SBRv055_34.16n';                                 % estm Navigation

%% Rinex to QM
WriteObs(obsfileBs);                                        % Base Station QMfile                 
renameBs = renameQMfile(obsfileBs);                    % Base Station's QMfile Renaming
WriteObs(obsfileRv);                                         % estm QMfile
renameRv = renameQMfile(obsfileRv);                    % estm's QMfile Renaming

%% �׹� RINEX ���� ����

%% ��Ÿ ���� : ����Ʈ ��ǥ ����
% TruePosBs = [];                                       % : Not exist True position
% TruePosRv = [];                                       % : Not exist True position
TruePosBs = [-3041235.578 4053941.677 3859881.013];   % : JPspace A point
TruePosRv = [-3041241.741 4053944.143 3859873.640];   % : JPspace B point
[YY, DOY] = obs2YYDOY(obsfileBs);

%% QM ���� �о�鿩�� ��ķ� �����ϰ�, ����� ����ġ ����
[arrQM1, FinalPRNs1, FinalTTs1] = ReadQM(renameBs);
QM1 = SelectQM(arrQM1, ObsType);
QM11 = SelectQM(arrQM1, ObsType2);
[arrQM2, FinalPRNs2, FinalTTs2] = ReadQM(renameRv);
QM2 = SelectQM(arrQM2, ObsType);
QM22 = SelectQM(arrQM2, ObsType2);

%% load QM
% load('DDevent.mat');

%% �׹��޽����� �о�鿩�� ��ķ� �����ϰ�, Klobuchar �� ����
eph = ReadEPH(navfileBs);
[al, be] = GetALBE(navfileBs);

%% ���̳ؽ� ���Ͽ��� �뷫���� ������ ��ǥ�� �̾Ƴ�
AppPos = GetAppPos(obsfileRv);
if AppPos(1) == 0
    AppPos = TruePosRv;
end
gd = xyz2gd(AppPos); AppLat = gd(1); AppLon = gd(2);

%% ���̳ؽ� ���Ͽ��� Base Station�� ��ǥ�� �����
Bs = PP(obsfileBs,navfileBs);              % without Correction
% Bs = PPwC(obsfileBs, navfileBs);            % with Correction
Bsxyz = [];
% [Bsgd, Bsgs, Bsutc, Bsla, Bslo, Bsh, Bsxyz] = GGA2gd(strcat(obsfileBs(1:(length(obsfileBs)-4)),'.ubx'));
% Bs = [Bsgs Bsxyz];

%% ���� �ð� sorting
% [year, month, days]= obs2date(obsfileBs);
% VST = [17, 12, 25];
% VET = [17, 14, 05];
% [gws, start_time] = date2gwgs(year, month, days, VST(1)-9, VST(2), VST(3)); start_time = round(start_time) - 17;
% [gws, end_time] = date2gwgs(year, month, days, VET(1)-9, VET(2), VET(3)); end_time = round(end_time) - 17;
% Bs = Bs(find(Bs(:,1) == start_time):find(Bs(:,1) == end_time),:);
%% True Distance
Truedis = 3.8;
%% �� QM ���Ͽ��� ����ð�(epoch) ����
if ~isempty(Bsxyz)
    FinalTTs = intersect(Bsgs(:, 1), QM2(:, 1));
else
    FinalTTs = intersect(Bs(:, 1), QM2(:, 1));
end

%% ������ �ʿ��� �ʱ�ġ ����
MaxIter = 4;
EpsStop = 1e-6;
x = AppPos';


%% �������� ����
NoEpochs = length(FinalTTs);
estm = zeros(NoEpochs, 6);  % : c1(gs), c2/3/4(X/Y/Z of Rv), c5(#Sats_obs), c6(#Sats_used)
nEst = 0;
cnt = 0;
No_Sat = 0;
for j = 1:NoEpochs
% for j = 28:28
    
    gs = FinalTTs(j);
    %% �ش� �ð� Bs�� ��ġ �� ã��
    Base(j,:) = Bs(find(Bs(:,1) == gs),1:4);
    TruePosBs_PP = Base(j,2:4);
    %% base�� gd
    base_gd(j,:) = xyz2gd(TruePosBs_PP(1:3)); % Base�� xyz�� gd �� ��ȯ
    AppLatBs = base_gd(j,1); AppLonBs = base_gd(j,2);
    %% �ش� �ð� gs�� ����ġ ���� �� ������� ���� ã��
    indexQM1 = find(QM1(:,1) == gs);
    QM1eBs = QM1(indexQM1,:);
    QM11eBs = QM11(indexQM1,:);
    %% rtklib�� ubx�� ��ȯ �� navfile�� ���������� �����Ҷ��� ����ϱ� ���� 
    existprnBs = intersect(unique(eph(:,18)), QM1eBs(:,2));
    arrSVBs = zeros(length(existprnBs),1);
    for kk = 1:length(existprnBs)
        arrSVBs(kk) = find(QM1eBs(:,2) == existprnBs(kk,1));
    end
    QM1eBs = QM1eBs(sort(arrSVBs),:);
    QM11eBs = QM11eBs(sort(arrSVBs),:);
    
    indexQM2 = find(QM2(:,1) == gs);
    QM1eRv = QM2(indexQM2,:);
    QM11eRv = QM22(indexQM2,:);
    %% rtklib�� ubx�� ��ȯ �� navfile�� ���������� �����Ҷ��� ����ϱ� ����
    existprnRv = intersect(unique(eph(:,18)), QM1eRv(:,2));
    arrSVRv = zeros(length(existprnRv),1);
    for kkk = 1:length(existprnRv)
        arrSVRv(kkk) = find(QM1eRv(:,2) == existprnRv(kkk,1));
    end
    QM1eRv = QM1eRv(sort(arrSVRv),:);
    QM11eRv = QM11eRv(sort(arrSVRv),:);
    
    
    Sats = intersect(QM1eBs(:, 2), QM1eRv(:, 2));
    NoSats = length(Sats); No_Sat = No_Sat + NoSats;
    
    %% �������� RS�� �ٸ����� OS ����/ SatsEl - c1(gs), c2(prn), c3(el)
    [SatsEl, indxRS] = PickRSel(gs, Sats, eph, TruePosBs_PP);  % : RS ��������
    RS = Sats(indxRS); RefSV(j,1) = RS;
    
    %% �������� ��ǥ ���� ��� - Bs ����
    icol = PickEPH(eph, RS ,gs);
    STT = GetSTTbrdc(gs, RS, eph, TruePosBs_PP);
    tc = gs -STT;
    vec_RS = GetSatPosNC(eph, icol, tc);
    vec_RS = RotSatPos(vec_RS, STT);
    S1BsRS = QM11eBs(find(QM11eBs(:,2) == RS), 4);
    S1RvRS = QM11eRv(find(QM11eRv(:,2) == RS), 4);
    for Iter = 1:MaxIter
        
        HTH = zeros(3,3);
        HTy = zeros(3,1);
        NoSatsUsed = NoSats;
        
        %% �� ������ ���� ����ġ, ���ġ, H��� ���
        for kS = 1:NoSats
            if kS == indxRS || SatsEl(kS, 3) < eleCut
                if SatsEl(kS, 3) < eleCut
                    NoSatsUsed = NoSatsUsed - 1;
                    disp([RS SatsEl(kS, 2)])
                end
                continue
            end
            cnt = cnt + 1;
            %% DD ����ġ ���� ��Ʈ --- ���Ŀ� for ���� ������ ���� �� 11/8/14
            OS = Sats(kS);
            S1BsOS = QM11eBs(find(QM11eBs(:,2) == OS), 4);
            S1RvOS = QM11eRv(find(QM11eRv(:,2) == OS), 4);
            obs_BsRS = QM1eBs(find(QM1eBs(:, 2) == RS), 4);
            obs_RvRS = QM1eRv(find(QM1eRv(:, 2) == RS), 4);
            obs_BsOS = QM1eBs(find(QM1eBs(:, 2) == OS), 4);
            obs_RvOS = QM1eRv(find(QM1eRv(:, 2) == OS), 4);
            obs = (obs_BsRS - obs_RvRS) - (obs_BsOS - obs_RvOS);
            %% DD ���ġ ���� ��Ʈ - ��Ÿ���� ��ǥ ���(�������� ��ǥ�� �̹� ��� �Ϸ�)
            icol = PickEPH(eph, OS, gs);
            STT = GetSTTbrdc(gs, OS, eph, x(1:3)'); % :  OS ������ġ�� ������ǥ �������� ���� 11/9/14
            tc = gs - STT;
            vec_OS = GetSatPosNC(eph, icol, tc);
            vec_OS = RotSatPos(vec_OS, STT);
            %% DD ���ġ ���� ��Ʈ - �Ÿ����ġ�� ���� ����� ���� DD ���ġ ���
            vec_BsRS = vec_RS - TruePosBs_PP;  com_BsRS = norm(vec_BsRS);
            vec_RvRS = vec_RS - x(1:3)';    com_RvRS = norm(vec_RvRS);
            vec_BsOS = vec_OS - TruePosBs_PP;  com_BsOS = norm(vec_BsOS);
            vec_RvOS = vec_OS - x(1:3)';    com_RvOS = norm(vec_RvOS);
            com = (com_BsRS - com_RvRS) - (com_BsOS - com_RvOS);
            y = obs -com;
            %% �� ���� az, el ����
            
            [azBsRS,elBsRS] = xyz2azel(vec_BsRS, AppLatBs, AppLonBs);
            azelBsRS(j,:) = [azBsRS, elBsRS];
            [azBsOS,elBsOS] = xyz2azel(vec_BsOS, AppLatBs, AppLonBs);
            azelBsOS(cnt,:) = [azBsOS, elBsOS];
            %% Weighting 
            S1RS(cnt,:) = [S1BsRS, S1RvRS,gs,OS];
            S1OS(cnt,:) = [S1BsOS, S1RvOS,gs,OS];
            DDel(cnt,:) = [elBsRS, elBsOS,gs,OS];
            W = DDMakeW_elsnr(S1RS(cnt,:),S1OS(cnt,:),DDel(cnt,:));
            weight(cnt,:) = W;
            W =1;
            %% H ��� ��� ��Ʈ
            H(1,1) = vec_RvRS(1)/com_RvRS - vec_RvOS(1)/com_RvOS;
            H(1,2) = vec_RvRS(2)/com_RvRS - vec_RvOS(2)/com_RvOS;
            H(1,3) = vec_RvRS(3)/com_RvRS - vec_RvOS(3)/com_RvOS;
            
            HTH = HTH + H'*W*H;
            HTy = HTy + H'*W*y;
        end
        
        xhat = inv(HTH) * HTy;
        x = x + xhat;
        
        if norm(xhat) < EpsStop;
            nEst = nEst + 1;
            estm(nEst,1) =gs;
            estm(nEst,2:4) =x;
            estm(nEst,5) = NoSats;
            estm(nEst,6) = NoSatsUsed;  % : ���������� �����ؾ� ��
            break;
        end
    end
    rover_gd(j,:) = xyz2gd(estm(j,2:4)); % rover�� xyz�� gd �� ��ȯ
    AppLat = rover_gd(j,1); AppLon = rover_gd(j,2);

end

%% �������� �м� & �׷��� �ۼ�
estm = estm(1:nEst, :);

% [dXYZ, dNEV] = PosTErrorsJOON(estm(:, 1), TruePosRv, estm(:, 2:5));         % A,B Point ���� ������
% [DDdXYZ, DDdNEV, DDdis, DDrms, DDstd] = PostErrorsDD(estm, Base, TruePosBs, TruePosRv); % A,B Point ���� ������
[DDdXYZ, DDdNEV, DDdis, DDrms, DDstd] = PostErrorsDDmv(estm, Base, Truedis,0,10);         % ������ ��ҿ��� �̵� ������  
% [QMnewBs, QMnewRv, QMnew] = DDSkyplot(QM1, QM2, eph, Base, estm);               % Skyplot
% DDPlotQM(renameBs, renameRv, 141)

%% Ư�� ���� �ð� Ž��
for TT = 1:length(estm(:,1))
    event = DDdis(TT,2);
    if event - Truedis > 1.5
        eT = estm(TT,1)-17;
        [yyy, mmo, ddd, hhh, mmm, sss] = gwgs2date(gws, eT);
        event_time(TT,:) = [hhh+9,mmm,sss,event-Truedis,1];
    else
        eT = estm(TT,1)-17;
        [yyy, mmo, ddd, hhh, mmm, sss] = gwgs2date(gws, eT);
        event_time(TT,:) = [hhh+9,mmm,sss,event-Truedis,0];
    end
end

% for i = 1:length(DDdis)
%     DDdis2D(i,:) = DDdis(i,1);
%     DDdis3D(i,:) = DDdis(i,2);
%     TrueDis(i,:) = DDdis(i,3);
%     figure(999)
%     hold on; grid on;
%     xlim([0 length(estm)])
%     plot(DDdis2D,'r-');
% %     plot(DDdis3D,'b-');
%     plot(TrueDis,'k-');
%     drawnow
% end
% 
% figure(999)
% hold on; grid on;
% xlabel({['Double Differencing dNE = ', num2str(DDrms(1)), '   std =', num2str(DDstd(1))]});
% % xlabel({['Double Differencing dNE = ', num2str(DDrms(1)), '   std =', num2str(DDstd(1))],...
% %     ['Double Differencing  3D = ', num2str(DDrms(2)), '   std =', num2str(DDstd(2))]});
% ylabel('Distance(meter)');
% legend('2D(dNE) Distance')
% % legend('2D(dNE) Distance','3D(d3D) Distance')