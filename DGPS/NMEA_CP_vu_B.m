close all; clear all;
tic
%% �Һ� ���� ����: ���� �ӵ�, ����ġ
CCC = 299792458.;   % CCC = Speed of Light [m/s]
%% �Ӱ���� ����
eleCut = 15;

%% NMEA ���Ϸ� ���� NMEAQM ����
POINT = 'B';                    % logginf Point      
DOY = '055'; YY = '16';         % Day of Year    
TIME = 'f';                     % Logging Hour(9a,10b,11c,12d,13e,14f,15g,16h,17i,18j,19k,
                                %             20l,21m,22n,23o,24p,1q,2r,3s,4t,5u,6v,7w,8x)  
DEVICE = 'vu';

NMEAfile = strcat('jpr',POINT,DOY,TIME,'_',DEVICE,'.txt');
NMEAQM = writeNMEA(NMEAfile);
% NMEAQM = writeNMEA2(NMEAfile);

%% PRC load
% PRCfile = 'JPRT160224.t1';
% PRCfile = 'JPRT160308.t1';
PRCfile = 'JPRT160224.t1';
[PRC_Sorted] = PRCsort(PRCfile, NMEAQM); 
% load('PRC_Sorted.mat');

%% NMEAQM �ڵ鸵
FinalTTs = intersect(unique(NMEAQM(:,1)), unique(PRC_Sorted(:,1))); 

%% �׹��޽����� �о�鿩�� ��ķ� �����ϰ�, Klobuchar �� ����, gps week, gps day ����
navfile = strcat('brdc',DOY,'0','.',YY,'n');
eph = ReadEPH(navfile);
[al, be] = GetALBE(navfile);
[gw, GD] = ydoy2gwgd(str2num(YY), str2num(DOY)); %: GPS WEEK ����

%% logging ������ ���� ������ ����
if POINT == 'A'
    TruePos = [-3041235.578 4053941.677 3859881.013];   % : JPspace A point
elseif POINT == 'B'
    TruePos = [-3041241.741 4053944.143 3859873.640];   % : JPspace B point
end

%% ������ �ʿ��� �ʱ�ġ ����
MaxIter = 4;
EpsStop = 1e-5;
ctr = 1;

%% �������� ����
NoEpochs = length(FinalTTs);
nEst = 0;
x_prc = zeros(4,1);

for j = 1:NoEpochs
    for Iter = 1:MaxIter
        HTH = zeros(4,4);
        HTy = zeros(4,1);
        HT_prc = zeros(4,1);      % pp-CP
                
        indexNMEAQM = find(NMEAQM(:,1) == FinalTTs(j));     % NMEA ��ü ��Ŀ��� ���� gs�� �� �� ����
        NMEAQM_1 = NMEAQM(indexNMEAQM,:);                   % NMEA ��ü ��Ŀ��� ����� gs ���� ������ ����
        indexPRC = find(PRC_Sorted(:,1) == FinalTTs(j));   % PRC ��ü ��Ŀ��� ���� gs�� �� �� ����
        PRC_1 = PRC_Sorted(indexPRC,:);                    % PRC ��ü ��Ŀ��� ����� gs ���� ������ ����
        
        NoSats = length(NMEAQM_1(:,5));                     % ���� gs�� ������
        gs = NMEAQM_1(1,1);                                 % ���� gs
        
        ddmm = [NMEAQM_1(1,2) NMEAQM_1(1,3) NMEAQM_1(1,4)]; % NMEA ��ü ��Ŀ��� ����� gs ���� Lat, Long, Alt ������ ����
        gd = ddmm2gd(ddmm);                                 % NMEA ��ü ��Ŀ��� ����� gs ���� ������ ����
        
        vec_site(j,:) = gd2xyz(gd);                              % Receiver xyz
        x = [vec_site(j,:) ctr]; x = x';
        ZHD = TropGPTh(vec_site(j,:), gw, gs);                 %: TROP: GPT
        
        for i = 1:NoSats
            prn = NMEAQM_1(i,5);                            % prn 
            prc = PRC_1(find(PRC_1(:,2) == prn),3);         % prn�� �´� PRC
            icol = PickEPH(eph, prn, gs);
            toe = eph(icol, 8); a = eph(icol, 19); b = eph(icol, 20); c = eph(icol, 21); Tgd = eph(icol, 23);
            
            %----- ��ȣ���޽ð� ���
            STT = GetSTTbrdc(gs, prn, eph, vec_site(j,:));
            tc = gs - STT;
                        
            %----- �����˵� ���
            vec_sat = GetSatPosNC(eph, icol, tc);          
            vec_sat = RotSatPos(vec_sat, STT);                %: �������� ���  
                                    
            %----- ���� RHO ���� ���
            vec_rho = vec_sat - vec_site(j,:);
            rho = norm(vec_rho);
            [az,el] = xyz2azel(vec_rho, gd(1), gd(2));        
                        
            if el >= eleCut %15
                W = 1;
%                 W = MakeW_elpr(el);
                
                dRel = GetRelBRDC(eph, icol, tc);
                dtSat = a + b*(tc - toe) + c*(tc - toe)^2 - Tgd + dRel;
                dIono = ionoKlob(al, be, gs, az, el, vec_site(j,:));
%                 dTrop_H = Hopfield(el, 0, vec_site(j,:), 9999);                   % Hopfield Model
%                 dTrop_S = Saastamoinen(el, 11, vec_site(j,:), 9999);               % Saastamoinen Model
                dTrop_G = ZHD2SHD(gw, gs, vec_site(j,:), el, ZHD);                   % GPT model
                y = CCC * (Tgd)+ dIono + dTrop_G ;
%                 cp = + dIono + dTrop_G + prc;                   % Correction-Projection
                cp =  + dIono + dTrop_G + prc;                   % Correction-Projection
%                 y =  + dIono + dTrop_G;
%                 cp = -CCC * (Tgd-dRel) - dIono - dTrop_G - prc;           
                H = [ -vec_rho(1)/rho -vec_rho(2)/rho -vec_rho(3)/rho 1];
                HTH = HTH + H'*W*H;
                HTy = HTy + H'*W*y;
                HT_prc = HT_prc + H'*W*cp;                                            % Correction-Projection
            end
                        
        end
%         xhat_prc =  inv(HTH) * HTy - inv(HTH) * HT_prc;                               % Correction-Projection
        xhat_prc = - inv(HTH) * HT_prc;                               % Correction-Projection
        x_prc = x + xhat_prc;                                         % Correction-Projection
                
%         if norm(xhat) < EpsStop;
        if Iter == 4;
            nEst = nEst + 1;
            estm(nEst,1) = gs;
            estm_prc(nEst,1) = gs;
            estm(nEst,2:5) = x(1:4);
            estm_prc(nEst,2:5) = x_prc(1:4);
            break;
        end

    end
    
    user_gd(j,:) = xyz2gd(estm_prc(j,2:4)); % user�� xyz�� gd �� ��ȯ
    user_xyz(j,:) = estm_prc(j,2:4);        % user xyz�� ��ķ� ��ȯ 
  
end

%% �������� �м� & �׷��� �ۼ�
estm = estm(1:nEst, :);

% [dXYZ, dNEV] = PosTErrorsJOON(estm(:, 1), TruePos, estm(:, 2:5));
% [prcdXYZ, prcdNEV] = PosTErrors2(estm_prc(:, 1), TruePos, estm_prc(:, 2:5));

[dXYZ, dNEV] = PosTErrors22(estm(:, 1), TruePos, estm(:, 2:5));
[dXYZ, dNEV] = PosTErrorsCP(estm(:, 1), TruePos, estm(:, 2:5), estm_prc(:, 2:5));
[prcdXYZ, prcdNEV] = PosTErrors2(estm_prc(:, 1), TruePos, estm_prc(:, 2:5));

standalone = [mean(sqrt(dNEV(:,1).^2+dNEV(:,2).^2)), mean(dNEV(:,3)), std(sqrt(dNEV(:,1).^2+dNEV(:,2).^2)), std(dNEV(:,3)),rms(sqrt(dNEV(:,1).^2+dNEV(:,2).^2)), rms(dNEV(:,3))]
cpcpcp = [mean(sqrt(prcdNEV(:,1).^2+prcdNEV(:,2).^2)), mean(prcdNEV(:,3)), std(sqrt(prcdNEV(:,1).^2+prcdNEV(:,2).^2)), std(prcdNEV(:,3)),rms(sqrt(prcdNEV(:,1).^2+prcdNEV(:,2).^2)), rms(prcdNEV(:,3))]

for mm = 1: length(estm(:,1))
    smartgd(mm,1:3) = xyz2gd(estm(mm,2:4));
    cpgd(mm,1:3) = xyz2gd(estm_prc(mm,2:4));
end

    