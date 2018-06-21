function [rename]= renameQMfile(filename)
%
%function [rename] = renameQMfile(obsfile)
%
%   Rename QMfile -> QplaceYYDOY
%   
%   input filename : Observation RINEX File
%
%   Example : renameQMfile('***.##o')
%
%   coded by Joonseong Gim, Jan 12, 2016
%   Modified by Joonseong Gim, Jan 15, 2016  % JPspace 표준 형식으로 renaming(ex:
%   Q'Place''YY''DOY'
%

fid = fopen(filename); % Obs file 
%fid = fopen('hyug118cd_.15o');
head_lines = 0;
while 1
   head_lines = head_lines+1;
   line = fgetl(fid);
   answer = findstr(line,'DATE');   % Obsfile에서 날짜가 있는 line 추출
   if ~isempty(answer), break;	end;
end;
Place = filename(1:4);
year = str2num(line(41:44));   
yy = str2num(line(43:44));
month = str2num(line(45:46));
day = str2num(line(47:48));
doy = date2doy(day,month,year); % year, month, day를 이용해 doy 계산
YY = num2str(yy); DOY = num2str(doy);
rename = strcat('Q', Place, YY, '0',DOY);
copyfile('QMfile', rename,'f');

