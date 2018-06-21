function [year, month, day, hour, minute, second]= obs2date(filename)
%
%function [rename] = renameQMfile(obsfile)
%
%   
%   input filename : Observation RINEX File
%
%   Example : renameQMfile('***.##o')
%
%   coded by Joonseong Gim, Jan 12, 2016
%   Modified by Joonseong Gim, Jan 15, 2016  % JPspace 표준 형식으로 renaming(ex:

%
fid = fopen(filename); % Obs file 
%fid = fopen('hyug118cd_.15o');
head_lines = 0;
while 1
   head_lines = head_lines+1;
   line = fgetl(fid);
   answer = findstr(line,'END OF HEADER');   % Obsfile에서 날짜가 있는 line 추출
   if ~isempty(answer), 
       head_lines = head_lines+1;
       line = fgetl(fid);
       break;	end;
end;
yr = str2num(line(2:3));   
if yr > 80
    year = yr + 1900;
else
    year = yr + 2000;
end

yy = str2num(line(2:3));
month = str2num(line(5:6));
day = str2num(line(8:9));
hour = str2num(line(11:12));
minute = str2num(line(14:15));
second = str2num(line(17:18));
doy = date2doy(day,month,year); % year, month, day를 이용해 doy 계산
YY = num2str(yy); DOY = num2str(doy);