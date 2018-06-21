function [estm] = PP_GPSGLO_SP3(SP3file, AppPos, QMfile)

%function [] = PP(obsfile, navfile)
%
%   Read the files(obs, nav) and Plot 'horizonal error', 'Vertical
%   error', from user_mean position(x,y,z)
%
%   input SP3file : Glonass SP3 file
%   input obsfile : Rinex Obsevation file
%   input QMfile : QMfile converted by T.I's converter
%
%   Example : PPwC_glo('QMfile')
%
%   coded by Joonseong Gim, Aug 18, 2016
%
% 
% clear all; close all;
% % 
% % %     load('QMjfr2');
% % QMfile = 'QMjfr2';
% QMfile = 'QSBSA_16050a';
% % obsfile = 'fr2r.obs';
% % 
% SP3file = 'igl18845.sp3';
% AppPos = [-3026795.499 4067267.161 3857084.459];
[arrSP3] = ReadSP3_GLO(SP3file);
%% �Һ� ���� ����: ���� �ӵ�, ����ġ
CCC = 299792458.;   % CCC = Speed of Light [m/s]
ObsType = 320;      % ����� ����ġ ���� - 120: C/A = C1
ObsType2 = 341;      % ����� ����ġ ���� - 141: snr = S1

real_x = -3041241.741;  % �뼺�������� ���� B ���� x
real_y = 4053944.143;   % �뼺�������� ���� B ���� y
real_z = 3859873.640;   % �뼺�������� ���� B ���� z
TruePos = [real_x real_y real_z];
%% �Ӱ���� ����
eleCut = 15;

%% rinex file �Է�
% obsfile = 'ubx33.obs';
% navfile = 'brdc0410.16n';
% obsfile = 'DAEB011c.16o';
% navfile = 'brdc0110.16n';

%% ������� QMfile�� obsevation Rinex ���� ���� ��¥, �ð����� ����
%     rename = renameQMfile(obsfile);
%     [YY, DOY] = obs2YYDOY(obsfile);
%     [gw, GD] = ydoy2gwgd(YY, DOY); %: GPS WEEK ����

%     [YY, DOY] = obs2YYDOY(obsfile);
%     ObsType = 120;
%     if DOY < 100
%         doy = strcat('0',num2str(DOY));
%         yy = num2str(YY);
%     else
%         doy = num2str(DOY);
%         yy = num2str(YY);
%     end

% navfile = strcat('brdc',doy,'0.',yy,'n');       % brdc file for GPS

%% QM ���� �о�鿩�� ��ķ� �����ϰ�, ����� ����ġ ����
[arrQM, FinalPRNs, FinalTTs] = ReadQM(QMfile);
QM = SelectQM(arrQM, ObsType);
QM2 = SelectQM(arrQM, ObsType2);
QM3 = SelectQM(arrQM, 131);
FinalTTs = FinalTTs(2:length(FinalTTs));
%% �׹��޽����� �о�鿩�� ��ķ� �����ϰ�, Klobuchar �� ����
%     eph = ReadEPH(navfile);
%     [al, be] = GetALBE(navfile);
%% ���̳ؽ� ���Ͽ��� �뷫���� ������ ��ǥ�� �̾Ƴ��� ���浵�� ��ȯ
% AppPos = GetAppPos(obsfile);
gd = xyz2gd(AppPos); AppLat = gd(1); AppLon = gd(2);

%% ������ �ʿ��� �ʱ�ġ ����
MaxIter = 4;
EpsStop = 1e-5;
ctr = 1;
x = [AppPos ctr]; x = x';

%% �������� ����
NoEpochs = length(FinalTTs);
nEst = 0;


for j = 1:NoEpochs
    indexQM = find(QM(:,1) == FinalTTs(j));
    QM_1 = QM(indexQM,:);
    QM_2 = QM2(indexQM,:);
    QM_3 = QM3(indexQM,:);
    
    for Iter = 1:MaxIter
        HTH = zeros(4,4);
        HTy = zeros(4,1);
        NoSats = length(QM_1);
        gs = QM_1(1,1);
        
        vec_site = x(1:3)';
        %                 ZHD = TropGPTh(vec_site, gw, gs);                 %: TROP: GPT
        
        for i = 1:NoSats
            prn = QM_1(i,2);
            obs = QM_1(i,4);
            S1 = QM_2(i,4);
            
          %% ���� ��ǥ ���
            STT = GetSTTsp3(gs, prn, arrSP3, x(1:3));
            tc = gs - STT;
            vec_sat = IntpSP3e1(arrSP3, prn, tc); vec_sat = vec_sat(2:4);
            %----- ���� RHO ���� ���
            vec_rho = vec_sat - vec_site;
            rho = norm(vec_rho);
            [az,el] = xyz2azel(vec_rho, AppLat, AppLon);
            if el >= eleCut %15
                %                         W = MakeW_elpr(el);
%                 W = MakeW_elsnr(el,S1);
                W = 1;
                %                         dRel = GetRelBRDC(eph, icol, tc);
                %                         dtSat = a + b*(tc - toe) + c*(tc - toe)^2 - Tgd + dRel;
                %                         dIono = ionoKlob(al, be, gs, az, el, vec_site);
                %                         dTrop_H = Hopfield(el, 11, vec_site, 9999);                   % Hopfield Model
                %                         dTrop_S = Saastamoinen(el, 11, vec_site, 9999);               % Saastamoinen Model
                %                         dTrop_G = ZHD2SHD(gw, gs, vec_site, el, ZHD);                   % GPT model
                %                         com = rho + x(4) - CCC * dtSat + dIono + dTrop_H;             % Hopfield Model
                %                         com = rho + x(4) - CCC * dtSat + dIono + dTrop_S;             % Saastamoinen Model
                %                         com = rho + x(4) - CCC * dtSat + dIono + dTrop_G;             % GPT Model
                %                         com = rho + x(4) - CCC * dtSat + dIono;                       % Satellite Retivistic
                com = rho + x(4);                               % Without Any correction
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
            estm(nEst,9) = NoSats;

            break;
        end
        
    end
    
    
end
%     end

estm = estm(find(estm(:,1) > 0),:);
for es = 1:length(estm(:,1))
    user_gd(es,:) = xyz2gd(estm(es,2:4)); % user�� xyz�� gd �� ��ȯ
    AppLat = user_gd(es,1); AppLon = user_gd(es,2);
    user_xyz(es,:) = estm(es,2:4);        % user xyz�� ��ķ� ��ȯ
end



