function [ObsSeq, NoTYPES] = GetObsSeq(obs_file)
% 
%function [ObsSeq, NoTYPES] = GetObsSeq(obs_file)
%
% DO: Extract the number of 'TYPES' and the sequence of them
%
% Copyright: Kwan-Dong Park@Koomin University
%
fid_obs = fopen(obs_file,'r');
ready = 0;
while ~ready
    s = fgetl(fid_obs);
    if length(s) > 72
%       TYPES Found        
        if s(61:72) == '# / TYPES OF'
            LineTYPES = s(1:length(s));
            ready = 1;
        end
%       TYPES Not Found        
        if s(61:72)=='END OF HEADE'
            ready = 1;
            fprintf('Error::: Not Found');
        end
    end
    
end
% NoTYPES = Number of TYPES of observables
NoTYPES = str2num(LineTYPES(6:6));

% DefTYPES = Default order of TYPES
% 1 = L1, 2 = L2, 3 = C1, 4 = P1, 5 = P2, 6 = D1, 7 = D2, 8 = S1, 9 = S2,
% 10 = C2
DefTYPES = char('L1','L2','C1','P1','P2','D1','D2','S1','S2','C2','C5','L5');

% Determine Observables Sequence
ObsSeq = zeros(1,NoTYPES);

for j = 1:NoTYPES
    Index1 = 11 + (j - 1) * 6;
    Index2 = Index1 + 1;
    T1 = LineTYPES(Index1:Index2);
    INDX_T1 = strmatch(T1,DefTYPES);
    ObsSeq(j) = INDX_T1; 
end

fclose(fid_obs);