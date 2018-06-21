%% �ڵ��ǻ�Ÿ� �������� �˰���
% 07/01/2016 : Joonseong

%% �Һ� ���� ���� : ���� �ӵ�, ����ġ
CCC = 299792458.;   % CCC = Speed of Light [m/s]
ObsType = 120;      % ����� ������ ���� - 120 : C/A = C1

%% �Ӱ���� ����
eleCut = 15;

%% QM ���� �ڵ鸵
FileQM1 = 'QZBLA_15334';
FileQM2 = 'QZBLB_15334';

%% �׹� RINEX ���� ����
FileNav = 'brdc3340.15n';

%% ��Ÿ ���� : ����Ʈ ��ǥ ���� & GPS Week
TruePosBs = [-3026789.236 4067255.523 3857098.106]; % : 15334 ����
TruePosRv = [-3026789.236 4067255.523 3857098.106]; % : 15334 ����

%% QM ���� �о�鿩�� ��ķ� �����ϰ�, ����� ����ġ ����
[arrQM1, FinalPRNs1, FinalTTs1] = ReadQM(FileQM1);
QM1 = SelectQM(arrQM1, ObsType);
[arrQM2, FinalPRNs2, FinalTTs2] = ReadQM(FileQM2);
QM2 = SelectQM(arrQM2, ObsType);

%% �� QM ���Ͽ��� ����ð�(epoch) ����
FinalTTs = intersect(QM1(:, 1), QM2(:, 1));

%% �׹��޽����� �о�鿩�� ��ķ� �����ϰ�, Klobuchar �� ����
eph = ReadEPH(FileNav);

%% ���̳ؽ� ���Ͽ��� �뷫���� ������ ��ǥ�� �̾Ƴ�
AppPos = TruePosRv; 

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
    %% �ش� �ð� gs�� ����ġ ���� �� ������� ���� ã��
    indexQM1 = find(QM1(:,1) == gs);
    QM1eBs = QM1(indexQM1,:);
    indexQM2 = find(QM2(:,1) == gs);
    QM1eRv = QM2(indexQM2,:);
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
[dXYZ, dNEV] = PosTErrorsJOON(estm(:, 1), TruePosRv, estm(:, 2:5));
% PosErrorsDD

            