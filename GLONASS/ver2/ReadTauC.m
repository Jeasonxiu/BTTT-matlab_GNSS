function [tau_c] = ReadTauC(GLON_file)

% function [tau_c] = ReadTauC(GLON_file)
%      Read TauC from glonass ephemerides file
% 
%      Coded by Miso Kim, October 14, 2014 

% GLON_file = 'brdc2760.14g';
fid = fopen(GLON_file, 'r');

ready = 0;

while ~ready
    s = fgetl(fid);
    if length(s) <= 80
        if s(61 : 64) == 'CORR'
            tau_c = str2num(s(23 : 40));
            break;
        end
    end
end