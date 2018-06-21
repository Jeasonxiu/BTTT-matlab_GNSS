function [year, month, day, hour, minute, second] = ReadRMC(s)
%
%function [year,month,day] = readRMC(string)
%
% DO: Extract Year-Month-Day from RMC
%
% <input>   string : RMC line
%
% <output>  Year(4-digit)/Month/Day : UTC
%
% <exception>
%           ReadRMC:NotEnoughData
%               수신 되다가 만 라인(파일 끝부분 발생)
%               <eg> $GPRMC,071509.0,A,3734.2632
%
% Copyright: Jinyi Kim, November 26, 2014@INHA University
%--- Modifications ---
% 11/27/2014 김진이가 만든 readRMC를 조금 변경함
%  1/12/2015 taeil Kim, <NotEnoughData> 예외 추가(1)
%  2/16/2015 taeil Kim, 시간도 추출하도록 변경

index=regexp(s(8:end),',','split'); % 콤마단위로 문자열 분리
if length(index)~=12                % 데이터 갯수 확인
    throw(MException('ReadRMC:NotEnoughData', ...
        'G#RMC:12개의 데이터중 %d개만 존재합니다.', length(index)));
end
%--- 날짜 추출 -------------------------------------------------------------
ddmmyy = index{9};
day   = str2num(ddmmyy(1:2));
month = str2num(ddmmyy(3:4));
year  = str2num(ddmmyy(5:end)) + 2000;
%--- 시간 추출 -------------------------------------------------------------
hhmmss = index{1};
hour  = str2num(hhmmss(1:2));
minute= str2num(hhmmss(3:4));
second= str2num(hhmmss(5:end));
end