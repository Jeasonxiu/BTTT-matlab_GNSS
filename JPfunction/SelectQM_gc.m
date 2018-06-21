function [SelectedQM] = SelectQM_gc(arrQM, g_ObsType, c_ObsType)
%
%function [SelectedQM] = SelectQM_gc(arrQM, g_ObsType, b_ObsType)
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
% 7/13/15 gps, bds ���ÿ� ó���ϵ��� ����
%
indxQM = find(arrQM(:,3) == g_ObsType);
indxQM2 = find(arrQM(:,3) == c_ObsType);
if length(indxQM) == 0 && length(indxQM2) == 0
    return
else
    g_QM = arrQM(indxQM,:);
    c_QM = arrQM(indxQM2,:);
end
g_QM(:,2)=g_QM(:,2)+100; % prn�� gps 100����, bds 200����� ����..
if c_ObsType > 200 & c_ObsType < 300
    c_QM(:,2)=c_QM(:,2)+200; % ���� WriteObs3 ����� �����ؾ�
elseif c_ObsType > 300
    c_QM(:,2)=c_QM(:,2)+300; % ���� WriteObs3 ����� �����ؾ�
end
SelectedQM=[g_QM; c_QM]; % gps, bds ������ ��ħ
SelectedQM=sortrows(SelectedQM,1); % �ð������� ����