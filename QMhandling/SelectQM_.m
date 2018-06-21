function [SelectedQM] = SelectQM(arrQM, ObsType)
%
%function [SelectedQM] = SelectQM(arrQM, ObsType)
%
% DO: Extract a specific ObsType from QM file
%
% <input>   arrQM
%           ObsType
%
% <output>  SelectedQM
%
% Copyright: Kwan-Dong Park @Jipyong Space
% -- Modifications --
% 9/27/14 코멘트(help디스플레이가 되지 않아 입출력 변수 설명 추가)
%
indxQM = find(arrQM(:,3) == ObsType);
if length(indxQM) == 0
    return
else
    SelectedQM = arrQM(indxQM,:);
end