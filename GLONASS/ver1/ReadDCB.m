function [DCB] = ReadDCB(ionex_file)

% IONEX 파일에서 DCB 읽어 DCB array로 저장함
% 2014.08.06 김미소

% GPS와 Glonass의 DCB를 모두 읽으며 gps = prn + 100 / glonass = prn + 200으로 구분.


% input 설정
% clc; clear all;
% ionex_file = 'igrg2120.14i';
% ionex_file = 'igsg2000.13i';
fid = fopen(ionex_file,'r');

ready = 0;
a = 0;

while ~ ready
    s = fgets(fid);
    if length(s) > 68
        if s(1:12) == 'DIFFERENTIAL' % 33번째행
            while ~ ready
                s = fgets(fid);
               
                if s(61:67) == 'STATION'
                    ready = 1;
                else
                    a = a + 1;
                    DCB(a,1) = str2num(s(5:6)); % GPS_PRN
                    DCB(a,2) = str2num(s(11:16)); % BIAS
                    DCB(a,3) = str2num(s(22:26)); % RMS
                    if s(4) == 'G'
                        DCB(a,1) = DCB(a,1) + 100;
                    else %s(4) == 'R'
                        DCB(a,1) = DCB(a,1) + 200;
                    end 
                end
                
            end
        end
        %     fclose(fid);
    end
end
