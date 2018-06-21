%% �ڵ��ǻ�Ÿ� �������� �˰���
% 07/01/2016 : Joonseong
clear all;

%% �Һ� ���� ���� : ���� �ӵ�, ����ġ
CCC = 299792458.;   % CCC = Speed of Light [m/s]
ObsType = 120;      % ����� ������ ���� - 120 : C/A = C1

%% �Ӱ���� ����
eleCut = 15;

%% QM ���� �ڵ鸵
obsfileBs = 'ihu1169a.11o';                                 % Reference Obsevation
navfileBs = 'brdc1690.11n';                                 % Reference Navigation
obsfileRv = 'ihu3169a.11o';                                 % estm Obsevation
navfileRv = 'brdc1690.11n';                                 % estm Navigation
WriteObs(obsfileBs);                                        % Base Station QMfile                 
renameBs = renameQMfile(obsfileBs);                    % Base Station's QMfile Renaming
WriteObs(obsfileRv);                                         % Rover QMfile
renameRv = renameQMfile(obsfileRv);                    % Rover's QMfile Renaming
%% �׹� RINEX ���� ����
% FileNav = 'brdc3340.15n';

%% ��Ÿ ���� : ����Ʈ ��ǥ ���� & GPS Week
% TruePosBs = [];                                       % : Not exist True position
% TruePosRv = [];                                       % : Not exist True position
% TruePosBs = [-3026789.236 4067255.523 3857098.106];   % : 15334 ����
% TruePosRv = [-3026789.236 4067255.523 3857098.106];   % : 15334 ����
TruePosBs = [-3026675.182 4067188.475 3857246.893];   % : ihu1 Phase Center
TruePosRv = [-3026675.978 4067187.900 3857246.933];   % : ihu3 Phase Center
% TruePosBs = [-3041235.578 4053941.677 3859881.013];   % : JPspace A point
% TruePosRv = [-3041241.741 4053944.143 3859873.640];   % : JPspace B point


%% QM ���� �о�鿩�� ��ķ� �����ϰ�, ����� ����ġ ����
[arrQM1, FinalPRNs1, FinalTTs1] = ReadQM(renameBs);
QM1 = SelectQM(arrQM1, ObsType);
[arrQM2, FinalPRNs2, FinalTTs2] = ReadQM(renameRv);
QM2 = SelectQM(arrQM2, ObsType);

%% �� QM ���Ͽ��� ����ð�(epoch) ����
FinalTTs = intersect(QM1(:, 1), QM2(:, 1));

%% �׹��޽����� �о�鿩�� ��ķ� �����ϰ�, Klobuchar �� ����
eph = ReadEPH(navfileBs);

%% ���̳ؽ� ���Ͽ��� �뷫���� ������ ��ǥ�� �̾Ƴ�
AppPos = TruePosRv; 

%% ���̳ؽ� ���Ͽ��� Base Station�� ��ǥ�� �����
Bs = PP(obsfileBs);

%% ������ �ʿ��� �ʱ�ġ ����
MaxIter = 4;
EpsStop = 1e-6;
x = AppPos';


%% �������� ����
NoEpochs = length(FinalTTs);
estm = zeros(NoEpochs, 6);  % : c1(gs), c2/3/4(X/Y/Z of Rv), c5(#Sats_obs), c6(#Sats_used)
nEst = 0;
for j = 1:NoEpochs
    
    gs = FinalTTs(j);
        %% �ش� �ð� Bs�� ��ġ �� ã��
    Base(j,:) = [gs,TruePosBs];
    %% �ش� �ð� gs�� ����ġ ���� �� ������� ���� ã��
    indexQM1 = find(QM1(:,1) == gs);
    QM1eBs = QM1(indexQM1,:);
    existprnBs = intersect(unique(eph(:,18)), QM1eBs(:,2));
    arrSVBs = zeros(length(existprnBs),1);
    for kk = 1:length(existprnBs)
        arrSVBs(kk) = find(QM1eBs(:,2) == existprnBs(kk,1));
    end
    QM1eBs = QM1eBs(sort(arrSVBs),:);
    indexQM2 = find(QM2(:,1) == gs);
    QM1eRv = QM2(indexQM2,:);
    existprnRv = intersect(unique(eph(:,18)), QM1eRv(:,2));
    arrSVRv = zeros(length(existprnRv),1);
    for kkk = 1:length(existprnRv)
        arrSVRv(kkk) = find(QM1eRv(:,2) == existprnRv(kkk,1));
    end
    QM1eRv = QM1eRv(sort(arrSVRv),:);
    Sats = intersect(QM1eBs(:, 2), QM1eRv(:, 2));
    NoSats = length(Sats);
    %% �������� RS�� �ٸ����� OS ����/ SatsEl - c1(gs), c2(prn), c3(el)
    [SatsEl, indxRS] = PickRSel(gs, Sats, eph, TruePosBs);  % : RS ��������
    RS = Sats(indxRS);
    %% �������� ��ǥ ���� ��� - Bs ����
    icol = PickEPH(eph, RS ,gs);
    STT = GetSTTbrdc(gs, RS, eph, TruePosBs);
    tc = gs -STT;
    vec_RS = GetSatPosNC(eph, icol, tc);
    vec_RS = RotSatPos(vec_RS, STT);
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
            %% DD ����ġ ���� ��Ʈ --- ���Ŀ� for ���� ������ ���� �� 11/8/14
            OS = Sats(kS);
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
            vec_BsRS = vec_RS - TruePosBs;  com_BsRS = norm(vec_BsRS);
            vec_RvRS = vec_RS - x(1:3)';    com_RvRS = norm(vec_RvRS);
            vec_BsOS = vec_OS - TruePosBs;  com_BsOS = norm(vec_BsOS);
            vec_RvOS = vec_OS - x(1:3)';    com_RvOS = norm(vec_RvOS);
            com = (com_BsRS - com_RvRS) - (com_BsOS - com_RvOS);
            y = obs -com;
            %% H ��� ��� ��Ʈ
            H(1,1) = vec_RvRS(1)/com_RvRS - vec_RvOS(1)/com_RvOS;
            H(1,2) = vec_RvRS(2)/com_RvRS - vec_RvOS(2)/com_RvOS;
            H(1,3) = vec_RvRS(3)/com_RvRS - vec_RvOS(3)/com_RvOS;
            
            HTH = HTH + H'*H;
            HTy = HTy + H'*y;
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
end

%% �������� �м� & �׷��� �ۼ�
estm = estm(1:nEst, :);
% [dXYZ, dNEV] = PosTErrorsJOON(estm(:, 1), TruePosRv, estm(:, 2:5));
% [DDdXYZ, DDdNEV, DDdis, DDrms, DDstd] = PostErrorsDD(estm, Base, TruePosBs, TruePosRv); 
[DDdXYZ, DDdNEV, DDdis, DDrms, DDstd] = PostErrorsDD3(estm, Bs, TruePosBs, TruePosRv); 

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