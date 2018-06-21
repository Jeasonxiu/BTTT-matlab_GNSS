clear; clc;
format long
% NMEA_file = 'PTCO4_joon.txt';
% NMEA1_file = 'PTCO4_joon_ori.txt';
% NMEA2_file = 'PTCO4_joon_adm.txt';
% NMEA_file = 'PTCO_CLT_170302.txt';
% NMEA1_file = 'PTCO_CLT_170302_ori.txt';
% NMEA2_file = 'PTCO_CLT_170302_adm.txt';
% NMEA_file = '160316.txt';
% NMEA1_file = 'SOND_170316_ori.txt';
% NMEA2_file = 'SOND_170316_adm.txt';
NMEA_file = 'DD1Rv.nmea';
NMEA1_file = 'DD1Rv_ori.txt';
NMEA2_file = 'DD1Rv_adm.txt';
% NMEA_file = 'DD_goGPS.txt';
% NMEA1_file = 'DD_goGPS_ori.txt';
% NMEA2_file = 'DD_goGPS_adm.txt';
yyyy = 2017; mo = 3; dd = 22; doy = 081;
fidNMEA = fopen(NMEA_file, 'r'); % ������ ����
fid_out1 = fopen(NMEA1_file, 'w'); % ������ ����
fid_out2 = fopen(NMEA2_file, 'w'); % ������ ����

line1=fgetl(fidNMEA); % ������ ù ���� ����
i = 0;
while(line1 ~= -1) % ���� ���� ������ �ݺ��� ����
    i = i + 1;
    line1=strtok(line1, '*'); % * ���ϸ� ���� (�ڵ��۾����� �ʿ���� �κ�)

    flag1=strread(line1(1:6),'%s'); % �ش� ���� �� 1~3��° ���� ����

    if strcmp(flag1, '$GPGGA')
        [h, m, s, x, y, z, la, lo, qi,nSats,ht] = NEWreadGGA2(line1);
        [gw, gs] = date2gwgs(yyyy, mo, dd, h, m, s);
        if ~isempty(x)
            NMEA(i, 1:10) = [h,m,s,la,lo,qi,gs,x,y,z];
            real_gs = gs + 18;
            fprintf(fid_out1,'%8.2f %13.10f %14.10f %6.3f %12.3f %12.3f %12.3f %1d \n',gs,la,lo,ht,x,y,z,qi);
            VRS_NMEA_adm(i,1) = gs;
            %         VRS_NMEA_adm(i,2:4) = cal_adm_pos([x,y,z],[doy,yyyy]);
            float_doy = yyyy + (doy-1)/365 + (h*3600 + m*60 + s)/(86400*365);
            VRS_NMEA_adm(i,2:4) = KPVeloTrans([x,y,z],float_doy);
            VRS_NMEA_gd(i,1) = gs;
            VRS_NMEA_gd(i,2:4) = xyz2gd(VRS_NMEA_adm(i,2:4));
            fprintf(fid_out2,'%8.2f %13.10f %14.10f %6.3f %12.3f %12.3f %12.3f %1d \n',gs,VRS_NMEA_gd(i,2),VRS_NMEA_gd(i,3),VRS_NMEA_gd(i,4),VRS_NMEA_adm(i,2),VRS_NMEA_adm(i,3),VRS_NMEA_adm(i,4),qi);
        end
    end
    line1=fgetl(fidNMEA); % ������ ù ���� ����
end
