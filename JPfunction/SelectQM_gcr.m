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
% 9/27/14 코멘트(help디스플레이가 되지 않아 입출력 변수 설명 추가)
% 7/13/15 gps, bds 동시에 처리하도록 수정, Hyunu Tae
% 8/24/15 gps, bds, glo 동시에 처리하도록 수정, Hyunu Tae

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
SelectedQM=[g_QM; c_QM; r_QM]; % gps, bds, glo 데이터 합침
SelectedQM=sortrows(SelectedQM,1); % 시간순으로 정렬