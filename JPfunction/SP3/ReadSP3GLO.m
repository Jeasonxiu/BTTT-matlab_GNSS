function [arrSP3] = ReadSP3GLO(file_sp3)
%
% ReadSP3_GLO('infile.sp3')
% GLONASS의 sp3파일을 읽어서 시간에 따라 행렬로 저장하는 코드
% result = [sec prn x(m) y(m) z(m) clk(10e-6s) yyyy mm dd hh min sec]
% 
% 원래의 ReadSP3를 GLONASS_SP3에 적용되도독 수정함
% 수정2 : 김미소
%        GPS의 SP3파일 -> GLONASS의 SP3를 읽도록 수정(ReadSP3 -> ReadSP3_GLO)


%% 함수 만들기 전 input 설정
% clc; clear all;
% file_sp3='igl18034.sp3';
% file_sp3 = 'Sta18034_sp3_glo.txt';
%
fid_sp3 = fopen(file_sp3,'r');
[s] = Get2ENDsp3(fid_sp3); %: 헤더 제거 -- 마지막줄 되돌리는 방법?
%% 행렬크기 미리 설정 위성최대갯수 24 X 96(15분이 24시간에 96개)
maxNsat = 24;
size_arrSP3 = maxNsat * 4;
arrSP3 = zeros(size_arrSP3, 6);
%% 루핑 시작
k = 0;
ready = 0;

while ~ready
    if s(1:3) == 'EOF'  %: 첫 열에 EOF가 기록된 경우 강제 종료한다.
        break;
    end
    
    if s(1,1) == '*'    %: Time-Tagging에는 "*"가 있음
        yr  = str2num(s(4:7));
        mon = str2num(s(9:10));
        day = str2num(s(12:13));
        hr  = str2num(s(15:16));
        min = str2num(s(18:19));
        sec = str2num(s(21:32));
        [gw, gs] = date2gwgs(yr, mon, day, hr, min, sec);
        
    elseif length(CutStr(s,' ',1)) == 2         %: PRN 1-9
        k = k + 1;
        prn = CutStr(s,' ',1);
        arrSP3(k,1) = round(gs);                %: time-tag (sec)
        arrSP3(k,2) = str2num(s(3:4));          %: PRN
        arrSP3(k,3) = str2num(s(6:18))*1000;    %: X (km->m)
        arrSP3(k,4) = str2num(s(19:32))*1000;   %: Y (km->m)
        arrSP3(k,5) = str2num(s(34:46))*1000;   %: Z (km->m)
        arrSP3(k,6) = str2num(s(47:60));        %: 시계오차 (micro-sec)
        
    else length(CutStr(s,' ',1)) == 4;          %: PRN 10-24
        k = k + 1;
        prn = CutStr(s,' ',1);
        arrSP3(k,1) = round(gs);                %: time-tag (sec)
        arrSP3(k,2) = str2num(s(3:4));          %: PRN
        arrSP3(k,3) = str2num(s(6:18))*1000;    %: X (km->m)
        arrSP3(k,4) = str2num(s(19:32))*1000;   %: Y (km->m)
        arrSP3(k,5) = str2num(s(34:46))*1000;   %: Z (km->m)
        arrSP3(k,6) = str2num(s(47:60));        %: 시계오차 (micro-sec)
        
    end
    
    s = fgets(fid_sp3);
    if feof(fid_sp3) == 1   %: 파일 끝에 도달하면 종료!
        ready = 1;
    end
end
fclose(fid_sp3);
