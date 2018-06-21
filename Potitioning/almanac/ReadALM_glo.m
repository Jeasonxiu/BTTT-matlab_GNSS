% function [gl_alm] = ReadALM_glo(alm_file)

%
% function [alm] = ReadALM_gl(alm_file)
%
% Read GLONASS ALM_FILE and return an array named 'gl_alm'
%
% July 04, 2008
%  
% [Hye-In Kim, Inha University]
alm_file ='alm_glo_17257.txt'; 

format long e;

%* Open Almanac File
fid_alm = fopen(alm_file,'r');

%* Define alm Array
gl_alm = zeros(25,16);

i = 0;
ready = 0;
while ready == 0
    s = fgets(fid_alm);
   if length(s) > 7
        i = i + 1;
        gl_alm(i,04) = str2num(s(1:3));     % day
        gl_alm(i,03) = str2num(s(4:5));     % month
        gl_alm(i,02) = str2num(s(7:10));   % year
        gl_alm(i,05) = str2num(s(13:18));  % the time of obtaining almanac from the 24 hours, with UTC
        s = fgets(fid_alm);
        gl_alm(i,01) = str2num(s(1:2))+300;    % PRN
        gl_alm(i,06) = str2num(s(5:6));    % frequency slot(-7-24)
        gl_alm(i,07) = str2num(s(9));       % health(0-1)
        gl_alm(i,08) = str2num(s(24:38));  % Equator time(sec)
        gl_alm(i,09) = str2num(s(41:55));  % correction GLONASS-UTC
        gl_alm(i,10) = str2num(s(58:72));  % correction GPS-GLONASS
        gl_alm(i,11) = str2num(s(75:89));  % correction of time KA GLONASS relative to the system time
        s = fgets(fid_alm);
        gl_alm(i,12) = str2num(s(1:14))*pi/180;  % the length of unit(radian)
        gl_alm(i,13) = str2num(s(17:29))*pi/180;  % 궤도 경사각 보정량(radian)
        gl_alm(i,14) = str2num(s(32:44))*pi/180;  % the argument of perigee(radian)
        gl_alm(i,15) = str2num(s(47:59));  % eccentricity
        gl_alm(i,16) = str2num(s(61:74));  % correction to the draconitic period
        gl_alm(i,17) = str2num(s(76:89));  % correction to the draconitic period   
        
        % Consider the SV Health
        if gl_alm(i,7) ~= 1 | gl_alm(i,01) == 6 | gl_alm(i,01) == 9 | gl_alm(i,01) == 24  
            gl_alm(i, :) = 0;
        end        
        
       else
        
        ready=1;
        
    end
end

gl_alm=gl_alm(1:i,:);

fclose(fid_alm);


