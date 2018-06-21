function [eph_glo] = ReadEPH_glo_kdpark(eph_file)
%
%function [eph_glo] = ReadEPH_glo(eph_file)
%
% DO: Read GLONASS navigation file (eg. brdc2120.14g) and array the content
%
% <input>   eph_file: GLONASS navigation file name
%
% <output>  eph_glo: Array of the ephemeris' content
%
%--- Modifications ---
% 12/30/14: 김미소 버전 정리
% 8/7/16: PRN에 300 더하기
% 12/25/16: c1(PRN), c2(gs)를 c1(gs), c2(PRN)로 변경
% 6/10/17: 인덱스 다시 정리. eGPS 연구노트 참고

%% 함수작성 전 입출력 테스트
% eph_file = 'brdc2120.14g';
%% 행렬 초기화
eph_glo = zeros(3000,19);
%% 항법 RINEX 파일 열기
fid_eph = fopen(eph_file, 'r');
%% END of Header 찾기
Get2END(fid_eph);
%%
i = 0;
ready = 0;
while ready == 0
    
    s = fgets(fid_eph);
    
    if length(s) > 7
        i = i + 1;
        eph_glo(i,2) = str2num(s(1:2)) + 300;   % PRN
        s([38 57 76]) = 'eee';
        
        yr  = str2num(s(4:5));
        mon = str2num(s(6:8));
        day = str2num(s(9:11));
        hr  = str2num(s(12:14));
        min = str2num(s(15:17));
        sec = str2num(s(18:20));
        
        if yr < 80
            yr = yr + 2000;
        end
        [gw, gs] = date2gwgs(yr, mon, day, hr, min, sec);
        
        eph_glo(i,01) = round(gs);          % time-sec
        eph_glo(i,12) = str2num(s(23:41));  % SV clk
        eph_glo(i,13) = str2num(s(42:60));  % SV relative frequency bias
        eph_glo(i,14) = str2num(s(61:79));  % message time frame
        
        s = fgets(fid_eph);
        s([19 38 57 76]) = 'eeee';
        eph_glo(i,03) = str2num(s(1:22));   % position X(km)
        eph_glo(i,06) = str2num(s(23:41));  % velocity X(km/s)
        eph_glo(i,09) = str2num(s(42:60));  % acceleration X(km/s^2)
        eph_glo(i,19) = str2num(s(61:79));  % health
        
        s = fgets(fid_eph);
        s([19 38 57 76]) = 'eeee';
        eph_glo(i,04) = str2num(s(1:22));   % position Y
        eph_glo(i,07) = str2num(s(23:41));  % velocity Y
        eph_glo(i,10) = str2num(s(42:60));  % acceleration Y
        eph_glo(i,16) = str2num(s(61:79));  % frequency number
        
        s = fgets(fid_eph);
        s([19 38 57 76]) = 'eeee';
        eph_glo(i,05) = str2num(s(1:22));   % position Z
        eph_glo(i,08) = str2num(s(23:41));  % velocity Z
        eph_glo(i,11) = str2num(s(42:60));  % acceleration Z
        eph_glo(i,17) = str2num(s(61:79));  % age of oper, information
        
    else
        ready=1;
    end
end
%%
eph_glo = eph_glo(1:i,:);
eph_glo(:,3:11) = eph_glo(:,3:11).*10^3; %: 단위변환--- km를 m로 변환
%% 깔끔한 뒷 정리
fclose(fid_eph);