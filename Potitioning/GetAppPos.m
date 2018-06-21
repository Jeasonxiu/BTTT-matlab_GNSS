function [AppPos] = GetAppPos(filename)
%
%function [AppPos] = GetAppPos(filename)
%
%   Read the given Observation RINEX file, make Approximation Position from file
%   
%   input filename : Observation RINEX File
%
%   Example : [AppPos] = GetAppPos('***.##o')
%
%   coded by Joonseong Gim, Jan 12, 2016
%
fid = fopen(filename);
head_lines = 0;
while 1
   head_lines = head_lines+1;
   line = fgetl(fid);
   answer = findstr(line,'APPROX POSITION XYZ');
   if ~isempty(answer), break;	end;
end;
App_x = str2double(line(1:14));
App_y = str2double(line(16:28));
App_z = str2double(line(30:42));
AppPos = [App_x App_y App_z];
