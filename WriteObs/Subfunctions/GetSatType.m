function [SatType] = GetSatType(obs_file)

% Copyright: Mi-So, Kim,  July 3, 2015

fid_obs = fopen(obs_file,'r');
ready = 0; n = 0;
while ~ready
    s = fgetl(fid_obs);
    if length(s) > 72
        if s(61:73) == 'SYS / # / OBS'
            if s(1) == 'G' || s(1) == 'C' || s(1) == 'R' 
            n = n + 1;
            SatType(n,1) = s(1);
            end
        end        
        if s(61:73) == 'END OF HEADER'
            ready = 1;
        end
    end % if length(s) > 72
end % while ~ ready
