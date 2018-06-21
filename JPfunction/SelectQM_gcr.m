function [SelectedQM] = SelectQM_gcr(arrQM, g_ObsType, c_ObsType, r_ObsType)
%
%function [SelectedQM] = SelectQM_gcr(arrQM, g_ObsType, c_ObsType, r_ObsType)
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
% 7/13/15 gps, bds ���ÿ� ó���ϵ��� ����, Hyunu Tae
% 8/24/15 gps, bds, glo ���ÿ� ó���ϵ��� ����, Hyunu Tae

indxQM = find(arrQM(:,3) == g_ObsType);
indxQM2 = find(arrQM(:,3) == c_ObsType);
indxQM3 = find(arrQM(:,3) == r_ObsType);
if length(indxQM) == 0 && length(indxQM2) == 0 && length(indxQM3) == 0
    return
else
    g_QM = arrQM(indxQM,:);
    c_QM = arrQM(indxQM2,:);
    r_QM = arrQM(indxQM3,:);
end
g_QM(:,2)=g_QM(:,2)+100; 
c_QM(:,2)=c_QM(:,2)+200; 
r_QM(:,2)=r_QM(:,2)+300; 
SelectedQM=[g_QM; c_QM; r_QM]; % gps, bds, glo ������ ��ħ
SelectedQM=sortrows(SelectedQM,1); % �ð������� ����