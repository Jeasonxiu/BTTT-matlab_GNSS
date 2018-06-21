function [eph] = ReadEPH(eph_file)
%
%function [eph] = ReadEPH(eph_file)
%
% DO:  Read ephemerides file and return an array named 'eph'
%
% Copyright: Kwan-Dong Park, February, 2001, Harvard-Smithsonian CfA
%--- Modifications ---
%  8/29/2001: To add a,b,c                  - now the size is 21
% 11/ 6/2011: To add SV health and Tgd      - now the size is 23
%  3/ 9/2014: To add IODE                   - now the size is 24
%

%% 파일 핸들과 결과 저장 행렬 초기화
fid_eph = fopen(eph_file,'r');
eph = zeros(3000,24);
%% 헤더 부분 패스
ready = 0;
while ~ready
    s = fgetl(fid_eph);
    if length(s) > 71
        if s(61:72) == 'END OF HEADE'
            ready = 1;
        end
    end
end
%% 한 줄씩 반복적으로 항법메시지 읽어 해당 칼럼에 저장
i = 0;
ready = 0;
while ready == 0
    s = fgets(fid_eph);
    if length(s) > 7
        i = i + 1;
        eph(i,18) = sscanf(s,'%d',1);   %PRN
        s([19 38 57 76]) = 'eeee';
        eph(i,19) = str2num(s(23:41));  %a
        eph(i,20) = str2num(s(42:60));  %b
        eph(i,21) = str2num(s(61:79));  %c

        s = fgets(fid_eph);
        s([19 38 57 76]) = 'eeee';
        eph(i,24) = str2num(s(1:22));   %IODE
        eph(i,01) = str2num(s(23:41));  %C_rs
        eph(i,02) = str2num(s(42:60));  %dn
        eph(i,03) = str2num(s(61:79));  %M_0

        s = fgets(fid_eph);
        s([19 38 57 76]) = 'eeee';
        eph(i,04) = str2num(s(1:22));   %C_uc
        eph(i,05) = str2num(s(23:41));  %e
        eph(i,06) = str2num(s(42:60));  %C_us
        eph(i,07) = str2num(s(61:79));  %sqtA

        s = fgets(fid_eph);
        s([19 38 57 76]) = 'eeee';
        eph(i,08) = str2num(s(1:22));   %t_oe
        eph(i,09) = str2num(s(23:41));  %C_ic
        eph(i,10) = str2num(s(42:60));  %Omega_0
        eph(i,11) = str2num(s(61:79));  %C_is

        s = fgets(fid_eph);
        s([19 38 57 76]) = 'eeee';
        eph(i,12) = str2num(s(1:22));   %i_0
        eph(i,13) = str2num(s(23:41));  %C_rc
        eph(i,14) = str2num(s(42:60));  %omega
        eph(i,15) = str2num(s(61:79));  %Omega_dot

        s = fgets(fid_eph);
        s([19 38 57 76]) = 'eeee';
        eph(i,16) = str2num(s(1:22));   %i_dot
        eph(i,17) = str2num(s(42:60));  %week_n
        
        s = fgets(fid_eph);
        s([19 38 57 76]) = 'eeee';
        eph(i,22) = str2num(s(23:41)) ; %SV health
        eph(i,23) = str2num(s(42:60)) ; %Tgd
        
        s = fgets(fid_eph);
    else        
        ready=1;
    end
end
%% 깔끔한 마무리
eph = eph(1:i,:);
eph(:,26) = eph(:,8) + 18;
fclose(fid_eph);
