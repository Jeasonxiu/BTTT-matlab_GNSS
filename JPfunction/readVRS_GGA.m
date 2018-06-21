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
fidNMEA = fopen(NMEA_file, 'r'); % 파일을 연다
fid_out1 = fopen(NMEA1_file, 'w'); % 파일을 연다
fid_out2 = fopen(NMEA2_file, 'w'); % 파일을 연다

line1=fgetl(fidNMEA); % 파일의 첫 줄을 읽음
i = 0;
while(line1 ~= -1) % 공백 줄을 읽으면 반복문 종료
    i = i + 1;
    line1=strtok(line1, '*'); % * 이하를 삭제 (코딩작업에서 필요없는 부분)

    flag1=strread(line1(1:6),'%s'); % 해당 줄의 앞 1~3번째 값을 읽음

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
    line1=fgetl(fidNMEA); % 파일의 첫 줄을 읽음
end
