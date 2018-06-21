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
% 9/27/14 �ڸ�Ʈ(help���÷��̰� ���� �ʾ� ����� ���� ���� �߰�)
%
indxQM = find(arrQM(:,3) == ObsType);
if length(indxQM) == 0
    return
else
    SelectedQM = arrQM(indxQM,:);
end