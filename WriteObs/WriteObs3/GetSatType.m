function [SatType] = GetSatType(obs_file)

% Copyright: Mi-So, Kim,  July 3, 2015
% --- Modifications, Hyunu, Tae, January 15, 2016
% 수정사항 : GAL/QZS/SBS 추가

% clc; clear all;
% obs_file = 'nnor2440.15o';

fid_obs = fopen(obs_file,'r');
ready = 0; n = 0;
while ~ready
    s = fgetl(fid_obs);
    if length(s) > 72
        if s(61:73) == 'SYS / # / OBS'
            if s(1) == 'G' || s(1) == 'C' || s(1) == 'R' || s(1) == 'J' || s(1) == 'E' || s(1) == 'S'
            n = n + 1;
            SatType(n,1) = s(1);
            end
        end        
        if s(61:73) == 'END OF HEADER'
            ready = 1;
        end
    end % if length(s) > 72
end % while ~ ready
