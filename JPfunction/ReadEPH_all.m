function [eph, getPRN, getSec] = ReadEPH_all(eph_file)
%
%function [eph] = ReadEPH_all(eph_file)
%
% DO:  Read ephemerides file and return an array named 'eph'
%
% Copyright: Kwan-Dong Park, February, 2001, Harvard-Smithsonian CfA
% Modification : Jinyi Kim, January, 2015, GPS Lab. INHA Univ.
%
%--- Modifications ---
%  8/29/2001: To add a,b,c                  - now the size is 21
% 11/ 6/2011: To add SV health and Tgd      - now the size is 23
%  3/ 9/2014: To add IODE                   - now the size is 24
% 14/ 1/2016: Galileo, QZSS, SBAS �߰� by ������
%

% < �ȳ�����>
% ���� GPS, GLO, BDS �� ���ؼ��� ���
% �߰� : Galileo, QZSS, SBAS ���

% < �������� >
% �� eph�� ���� 24������ 25���� �����ȴ�.
%    BDS�� GALILEO ��� TGD�� 2���̱� ������ 23, 24���� ����
%    ���� IODE�� 24�� -> 25�� �� ����

% �� QZSS    : GPS�� ����
%    GLONASS,SBAS : GPS�� ü�谡 ���� �ٸ�!!!
%    GALILEO : TGD �ٸ�, Week Number�� ������ �ٸ��� ��ġ�� ����
%    BEIDOU  : TGD �ٸ�, Week Number�� ����� ��ġ�� �ٸ� (GPS = BDT + 1356)
%    ex) GPS WN = 1826, QZS WN = 1826, GAL WN = 1826, BDT WN = 470

% clear all;
% % eph_file = 'brdm2750.15p';
% eph_file = 'brdm0250_.17p';
% eph_file = 'brdm0100.15p';
%% ���� �ڵ�� ��� ���� ��� �ʱ�ȭ
fid_eph = fopen(eph_file,'r');
eph = zeros(3000,26);
%% ��� �κ� �߱��~
ready = 0;
while ~ready
    s = fgetl(fid_eph);
    if length(s) > 71
        if s(61:72) == 'END OF HEADE'
            ready = 1;
        end
    end
end
%% �� �پ� �ݺ������� �׹��޽��� �о� �ش� Į���� ����
i = 0;
ready = 0;
getPRN=0; getSec=0;
while ready == 0
    s = fgets(fid_eph);
    if length(s) > 7                    %[G,E,C,J]  %[R,S]
        i = i + 1;
        SSI = sscanf(s, '%c',1);
        NN = str2num(s(2:3));
        if ~isempty(NN) 
            eph(i,18) = PRNnumber(SSI,NN);  %PRN        %
        end
        %%
        yr  = str2num(s(7:8));
        mon = str2num(s(9:11));
        day = str2num(s(12:14));
        hr  = str2num(s(15:17));
        min = str2num(s(18:20));
        sec = str2num(s(21:23));
        if yr < 50
            yr = yr + 2000;
        end
        [gw, gs] = date2gwgs(yr, mon, day, hr, min, sec);
        eph(i,26) = round(gs)+LeapSeconds([yr, mon, day, hr, min, sec]);  %Toc �߰� 
        %%
        
        s([20 39 58 77]) = 'eeee';
        if SSI == 'R'|| SSI == 'S'
            eph(i,1) = PRNnumber(SSI,NN);  %PRN        %
            eph(i,12) = str2num(s(24:42));  % SV clk
            eph(i,13) = str2num(s(43:61));  % SV relative frequency bias
            eph(i,02) = round(gs);   % t_oe
            
            s = fgets(fid_eph);
            s([20 39 58 77]) = 'eeee';
            eph(i,03) = str2num(s( 1:23));  % position X(km)
            eph(i,06) = str2num(s(24:42));  % velocity X(km/s)
            eph(i,09) = str2num(s(43:61));  % acceleration X(km/s^2)
            eph(i,15) = str2num(s(62:80));  % SV health
            
            s = fgets(fid_eph);
            s([20 39 58 77]) = 'eeee';
            eph(i,04) = str2num(s( 1:23));  % position Y
            eph(i,07) = str2num(s(24:42));  % velocity Y
            eph(i,10) = str2num(s(43:61));  % acceleration Y
            eph(i,16) = str2num(s(62:80));  % frequency number
            
            s = fgets(fid_eph);
            s([20 39 58 77]) = 'eeee';
            eph(i,05) = str2num(s( 1:23));  % position Z
            eph(i,08) = str2num(s(24:42));  % velocity Z
            eph(i,11) = str2num(s(43:61));  % accerleration Z
            eph(i,14) = str2num(s(62:80));  % age if oper, information
            eph(i,3:8) = eph(i,3:8).*10^3;
            eph(i,9:11) = eph(i,9:11).*10^3;
        elseif SSI == 'I'
            eph(i,19) = str2num(s(24:42));  %a          %
            eph(i,20) = str2num(s(43:61));  %b          %
            eph(i,21) = str2num(s(62:80));  %c          %
            
            s = fgets(fid_eph);
            s([20 39 58 77]) = 'eeee';
            eph(i,25) = str2num(s( 1:23));  %IODE       %
            eph(i,01) = str2num(s(24:42));  %C_rs       %
            eph(i,02) = str2num(s(43:61));  %dn         %
            eph(i,03) = str2num(s(62:80));  %M_0        %
            
            s = fgets(fid_eph);
            s([20 39 58 77]) = 'eeee';
            eph(i,04) = str2num(s( 1:23));  %C_uc       %
            eph(i,05) = str2num(s(24:42));  %e          %
            eph(i,06) = str2num(s(43:61));  %C_us       %
            eph(i,07) = str2num(s(62:80));  %sqtA       %
            
            s = fgets(fid_eph);
            s([20 39 58 77]) = 'eeee';
            eph(i,08) = str2num(s( 1:23));  %t_oe       %
            eph(i,09) = str2num(s(24:42));  %C_ic       %
            eph(i,10) = str2num(s(43:61));  %Omega_0    %
            eph(i,11) = str2num(s(62:80));  %C_is       %
            
            s = fgets(fid_eph);
            s([20 39 58 77]) = 'eeee';
            eph(i,12) = str2num(s( 1:23));  %i_0
            eph(i,13) = str2num(s(24:42));  %C_rc
            eph(i,14) = str2num(s(43:61));  %omega
            eph(i,15) = str2num(s(62:80));  %Omega_dot
            
            s = fgets(fid_eph);
            s([20 39 58 77]) = 'eeee';
            eph(i,16) = str2num(s( 1:23));  %i_dot
            eph(i,17) = str2num(s(43:61));  %week_n
            
            s = fgets(fid_eph);
            s([20 39 58 77]) = 'eeee';
            eph(i,22) = str2num(s(24:42)) ; %SV health
            eph(i,23) = str2num(s(43:61)) ; %[G,J] Tgd  %[C] TGD1 B1/B3
            %[E] BGD E5a/E1
            eph(i,24) = 0 ; %[G,J] IODC %[C] TGD1 B1/B3
            %[E] BGD E5a/E1
            s = fgets(fid_eph);
        else
            eph(i,19) = str2num(s(24:42));  %a          %
            eph(i,20) = str2num(s(43:61));  %b          %
            eph(i,21) = str2num(s(62:80));  %c          %
            
            s = fgets(fid_eph);
            s([20 39 58 77]) = 'eeee';
            eph(i,25) = str2num(s( 1:23));  %IODE       %
            eph(i,01) = str2num(s(24:42));  %C_rs       %
            eph(i,02) = str2num(s(43:61));  %dn         %
            eph(i,03) = str2num(s(62:80));  %M_0        %
            
            s = fgets(fid_eph);
            s([20 39 58 77]) = 'eeee';
            eph(i,04) = str2num(s( 1:23));  %C_uc       %
            eph(i,05) = str2num(s(24:42));  %e          %
            eph(i,06) = str2num(s(43:61));  %C_us       %
            eph(i,07) = str2num(s(62:80));  %sqtA       %
            
            s = fgets(fid_eph);
            s([20 39 58 77]) = 'eeee';
            eph(i,08) = str2num(s( 1:23));  %t_oe       %
            eph(i,09) = str2num(s(24:42));  %C_ic       %
            eph(i,10) = str2num(s(43:61));  %Omega_0    %
            eph(i,11) = str2num(s(62:80));  %C_is       %
            
            s = fgets(fid_eph);
            s([20 39 58 77]) = 'eeee';
            eph(i,12) = str2num(s( 1:23));  %i_0
            eph(i,13) = str2num(s(24:42));  %C_rc
            eph(i,14) = str2num(s(43:61));  %omega
            eph(i,15) = str2num(s(62:80));  %Omega_dot
            
            s = fgets(fid_eph);
            s([20 39 58 77]) = 'eeee';
            eph(i,16) = str2num(s( 1:23));  %i_dot
            eph(i,17) = str2num(s(43:61));  %week_n
            
            s = fgets(fid_eph);
            s([20 39 58 77]) = 'eeee';
            eph(i,22) = str2num(s(24:42)) ; %SV health
            eph(i,23) = str2num(s(43:61)) ; %[G,J] Tgd  %[C] TGD1 B1/B3
            %[E] BGD E5a/E1
            eph(i,24) = str2num(s(62:80)) ; %[G,J] IODC %[C] TGD1 B1/B3
            %[E] BGD E5a/E1
            
            
            
            %% ����ó��... ���� health ���
            if eph(i,22) ~= 0 && eph(i,18)<400
                getSec(i)=eph(i,8);
                eph(i,8) = 0;
                getPRN(i)=eph(i,18);
            end
            %%
            s = fgets(fid_eph);
        end
        
    else
        ready=1;
    end
end
%% ����ó��... �˵��� ������ �̻�� ó��
for i=1:length(eph)-1
    if eph(i,3)==eph(i+1,3)
        eph(i+1,8)=0;
    end
end
%% GPS, GLONASS, BEIDOU�� �׹� �޽����� ��� (���� ���� ����)
GPS = PRNnumber('G',0)/100; indexGPS = find(floor(eph(:,18)/100) == GPS);
GLO = PRNnumber('R',0)/100; indexGLO = find(floor(eph(:,18)/100) == GLO);
BDS = PRNnumber('C',0)/100; indexBDS = find(floor(eph(:,18)/100) == BDS);
QZS = PRNnumber('J',0)/100; indexQZS = find(floor(eph(:,18)/100) == QZS);
GAL = PRNnumber('E',0)/100; indexGAL = find(floor(eph(:,18)/100) == GAL);
SBS = PRNnumber('S',0)/100; indexSBS = find(floor(eph(:,18)/100) == SBS);
indexGNS = union(indexGPS,indexBDS); indexGNS = union(indexGNS,indexGLO);
indexGNS = union(indexGNS,indexGAL); indexGNS = union(indexGNS,indexQZS);
indexGNS = union(indexGNS,indexSBS);
eph = eph(indexGNS,:); i = length(eph);

%% ����� ������
% eph = eph(1:i,:);


getPRN=unique(getPRN);
getSec=unique(getSec);
fclose(fid_eph);
