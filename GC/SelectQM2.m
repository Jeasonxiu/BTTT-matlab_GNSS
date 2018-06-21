function [out] = SelectQM2(arrQM, prn, obs1, obs2)
%
% function [out] = SelectQM2(arrQM, prn, obs1, obs2)
% 
% DO: Extracts two types of observalble from the input QM array
%
% <input>   - arrQM: Quick Measurements "array"
%           - prn: Satellite PRN ID
%           - obs1/obs2: Two types of observables to extract
%
% <output>  - 'out' array with c1(Time), c2(PRN), c3(obs1), c4(obs2)
%
% Copyright: Kwan-Dong Park, 13-10-01@LDEO
%---- Modifications
% 2/21/14: c1에 시각, c2에 PRN을 저장하도록 변경함. 이전에는 c1(PRN), c2(TT)

%%
colPRN = 2;
colObs = 3;

%-- Assign an empty array to 'out'
out = [];

%-- Find Satellite
indexPRN = find(arrQM(:,colPRN) == prn);
if isempty(indexPRN)
    fprintf('\n** ERROR ** Satellite %2d : NOT found\n\n', prn)
    return;
end
arrPRN = arrQM(indexPRN,:);
TTs = unique(arrPRN(:,1));
nEpochs = length(TTs);
out = zeros(nEpochs, 4);        
nOut = 0;
for k = 1:nEpochs
    indexTT = find(arrPRN(:,1) == TTs(k));
    arrPRNTT = arrPRN(indexTT,:);
    indexObs1 = find(arrPRNTT(:,colObs) == obs1);
    indexObs2 = find(arrPRNTT(:,colObs) == obs2);
    if isempty(indexObs1) || isempty(indexObs2)
        fprintf('** Warning: either %4d or %4d missing\n', obs1, obs2)
        continue;
    end
    nOut = nOut + 1;
    o1 = arrPRNTT(indexObs1, 4);
    o2 = arrPRNTT(indexObs2, 4);
    out(nOut,1) = TTs(k);
    out(nOut,2) = prn;
    out(nOut,3) = o1;
    out(nOut,4) = o2;
    %-- fprintf('%2d %10.2f %15.3f %15.3f\n', TTs(k), prn, o1, o2) 
end
out = out(1:nOut,:);
