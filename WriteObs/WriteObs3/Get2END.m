function [] = Get2END(fid_obs)
%
%function [] = Get2END(fid_obs)
% 
% Get to the last part of observation RINEX Header: 'END OF HEADER'
%
% <input> fid_obs: Input file FILE ID(handler)
%
% Copyright: Kwan-Dong Park, July 22, 2006@Inha University
%
ready = 0;
while ready == 0
    s = fgets(fid_obs);
    if length(s) > 71
        if s(61:71) == 'END OF HEAD'
            ready = 1;
        end
    end
end
