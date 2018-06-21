function [eph_glo] = ReadEPH_glo(eph_file)

%function [eph] = ReadEPH(eph_file)
%
%   Read ephemerides file and return an array named 'eph'
%
%       - Coded by Kwan-Dong Park, February, 2001, Harvard-Smithsonian CfA


% clc; clear all;

% 함수만들기 전 input 설정
% eph_file = 'MCCT1210_14G.txt';

%* Open Epheremides File
fid_eph = fopen(eph_file,'r');

%* END of Header
Get2END(fid_eph);

%* Define EPH Array
eph_glo = zeros(3000,17);

i = 0;
ready = 0;
while ready == 0
    
    s = fgets(fid_eph);
    
    if length(s) > 7
        i = i + 1;
        eph_glo(i,1) = str2num(s(1:2));   % PRN
        s([38 57 76]) = 'eee';
        
        yr  = str2num(s(4:5));
        mon = str2num(s(6:8));
        day = str2num(s(9:11));
        hr  = str2num(s(12:14));
        min = str2num(s(15:17));
        sec = str2num(s(18:20));
        if yr < 50
            yr = yr + 2000;
        end
        [gw, gs] = date2gwgs(yr, mon, day, hr, min, sec);
        
        eph_glo(i,02) = round(gs);   % time-sec
        eph_glo(i,12) = str2num(s(23:41));  % SV clk
        eph_glo(i,13) = str2num(s(42:60));  % SV relative frequency bias
        eph_glo(i,14) = str2num(s(61:79));  % massage time frame
        
        s = fgets(fid_eph);
        s([19 38 57 76]) = 'eeee';
        eph_glo(i,03) = str2num(s(1:22));   % position X(km)
        eph_glo(i,06) = str2num(s(23:41));  % velocity X(km/s)
        eph_glo(i,09) = str2num(s(42:60));  % acceleration X(km/s^2)
        eph_glo(i,15) = str2num(s(61:79));  % Health
        
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

eph_glo=eph_glo(1:i,:);
eph_glo(:,3:11) = eph_glo(:,3:11).*10^3; % km -> m
 
fclose(fid_eph);