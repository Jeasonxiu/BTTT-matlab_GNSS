function [SelectedQM] = SelectQM(arrQM, ObsType)
%
indexQM = find(arrQM(:,3) == ObsType);
if length(indexQM) == 0
    return
else
    SelectedQM = arrQM(indexQM,:);
end