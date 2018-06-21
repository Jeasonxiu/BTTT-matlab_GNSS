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
% 9/27/14 코멘트(help디스플레이가 되지 않아 입출력 변수 설명 추가)
% 7/13/15 gps, bds 동시에 처리하도록 수정
%
indxQM = find(arrQM(:,3) == g_ObsType);
indxQM2 = find(arrQM(:,3) == c_ObsType);
if length(indxQM) == 0 && length(indxQM2) == 0
    return
else
    g_QM = arrQM(indxQM,:);
    c_QM = arrQM(indxQM2,:);
end
g_QM(:,2)=g_QM(:,2)+100; % prn을 gps 100번대, bds 200번대로 수정..
if c_ObsType > 200 & c_ObsType < 300
    c_QM(:,2)=c_QM(:,2)+200; % 차후 WriteObs3 변경시 수정해야
elseif c_ObsType > 300
    c_QM(:,2)=c_QM(:,2)+300; % 차후 WriteObs3 변경시 수정해야
end
SelectedQM=[g_QM; c_QM]; % gps, bds 데이터 합침
SelectedQM=sortrows(SelectedQM,1); % 시간순으로 정렬