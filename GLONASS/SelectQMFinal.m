function [SelectedQM, FinalPRNs, FinalTTs] = SelectQM(arrQM, ObsType)
%
%function [SelectedQM] = SelectQM(arrQM, ObsType)
%
% DO: Extract a specific ObsType from QM file plus Final TT and PRNS
%
% <input>   arrQM
%           ObsType
%
% <output>  SelectedQM
%           FinalPRNs
%           FinalTTs
%
% Copyright: Kwan-Dong Park; June 10, 2017, @PPSoln
% -- Modifications --
% 
%
indxQM = find(arrQM(:,3) == ObsType);
if length(indxQM) == 0
    return
else
    SelectedQM = arrQM(indxQM,:);
end
%% 
FinalPRNs = unique(SelectedQM(:,2));
FinalTTs  = unique(SelectedQM(:,1));