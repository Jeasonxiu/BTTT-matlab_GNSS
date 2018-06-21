function [DeletedQM, FinalPRNs, FinalTTs] = DelQMprn(arrQM, prns)
%
%function [DeletedQM] = DelQMprn(arrQM, prns)
%
% DO: Delete a specified set of PRNs from QM file
%
% <input>   arrQM: Input QM file; c1(TT), c2(PRN), c3(ObsType), c4(Obs)
%           prns: a vecotor PRN IDs
%
% <output>  DeletedQM
%
% Copyright: Kwan-Dong Park @Jipyong Space, September 30, 2014
% -- Modifications --
% ?/?/?
%
%% 삭제할 PRN이 입력된 QM에 존재하지 않을 경우에 대비
DeletedQM = arrQM;
%% 
for k = 1:length(prns)
    prn = prns(k);
    indxQM = find(arrQM(:,2) ~= prn);
    if length(indxQM) == length(arrQM) 
        fprintf(1,'** Warning: PRN %2d not found\n', prn)
    else
        DeletedQM = arrQM(indxQM,:);
    end
    arrQM = DeletedQM;
end
%% ReadQM과 같은 OUTPUT 출력
PRNs = DeletedQM(:,2);
TTs = DeletedQM(:,1);
FinalPRNs = unique(PRNs);
FinalTTs = unique(TTs);
