function [OL] = csCombPL(qm_file, prn, obs_PR, obs_CP)
%
%function [OL] = csCombPL(qm_file, prn, obs_PR, obs_CP)
%
% DO: Combine Pseudo-Range and Phase measurements for Cycle Slip detection
% 
% <input>   QM_file: QM File 
%           prn: Satellite PRN 
%           obs_PR: Pseudo-range type   PR: 120(C1), 121(P1), 122(P2)
%           obs_CR: Carrier-Phase type  CP: 111(L1), 112(L2)
% 
% <output>  OL: Outliers    c1(gs); c2(PRN); c3(Delta N)
%           Figure
%
% Copyright: Kwan-Dong Park, 13-09-27 @LDEO
% ----- Modificattion
%  2/21/14-- c1/c2 순서 변경 - SelectQM2 참고
% 11/17/14-- C1(P1)&L1 이외에도 L2C(P2)&L2 조합도 CS 검출하도록 수정

%% CP의 종류에 따라 LAMBDA 변경하도록 수정 11/17/14
CCC = 299792458.;
if obs_CP == 111
    freq = 1575.42e6;
elseif obs_CP == 112
    freq = 1227.60e6;
else
    Warning('Check the observation type of CP : 111 or 112')
end
lamb = CCC/freq;

%% Format of CS: c1(TT), c2(PRN), c3(PR), c4(PH)
colT = 1;
colPR = 3;
colPH = 4;
%%
[arrQM, dummy, dummy] = ReadQM(qm_file);
cs = SelectQM2(arrQM, prn, obs_PR, obs_CP);
if isempty(cs)
    OL = [];
    return;
end
%% Set Interval to the most frequent occurence of TD
tintv = mode(diff(cs(:,colT)));
%% PR/PH 차이값 계산
nEpochs = length(cs);
comb = zeros(nEpochs,3);
nComb = 0;
for k = 2:nEpochs
    td   = cs(k, colT) - cs(k-1, colT);
    if td == tintv
        nComb = nComb + 1;
        td_PR = cs(k, colPR) - cs(k-1, colPR);  %td: time difference
        td_PH = cs(k, colPH) - cs(k-1, colPH);
        td_N = td_PH - td_PR/lamb;
        % fprintf('%3d %15.3f %15.3f %15.3f \n', td, td_P, td_L, td_N)
        comb(nComb, 1) = cs(k, colT);
        comb(nComb, 2) = prn;
        comb(nComb, 3) = td_N;
    end
end
%% Format of COMB: c1(PRN), c2(TT), c3(COMB)
comb = comb(1:nComb,:);
%% Find OUTLIERS
limitOL = 20;
indexOL = find(abs(comb(:,3)) > limitOL);    %: OL: outliers
OL = comb(indexOL,:);
%* Plot
tHour = gs2h24(comb(:,colT));
yObs = comb(:,3);
plot(tHour, yObs, 'o:'); axis([0 24 min(yObs) max(yObs)]);
xlabel('Hour'), ylabel('Code-Phase Comb (m)'), title(num2str(prn))
