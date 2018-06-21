function [] = csCombPL_all(qm_file, obs_PR, obs_CP)

%function [] = csCombPL_all(qm_file, obs_PR, obs_CP)
%
% DO: Multiple Code-Carrier combinations for Cycle Slip detection
% 
% <input>   QM_file 
%           obs_PR: Pseudo-range type   PR: 120(C1), 121(P1), 122(P2)
%           obs_CR: Carrier-Phase type  CP: 111(L1), 112(L2)
% 
% output: Figure
%
% Copyright: Kwan-Dong Park, 11/19/2014 @Jipyong Space

%% QM 파일 읽기
[arrQM, FinalPRNs, FinalTTs] = ReadQM(qm_file);
FinalPRNs = unique(arrQM(find(arrQM(:,3) < 200),2));
tS = gs2h24(FinalTTs(1)); tE = gs2h24(FinalTTs(end));
yMax = 20; 
%% csMW_one으로 각 위성에 대해 분석 - 그래프 그리기
nFigs_page = 4;
for kS = 1:length(FinalPRNs)
    prn = FinalPRNs(kS);
    %% subplot과 figure 번호 산출
    [nPage, subp] = DivideSubplot(nFigs_page, kS);
    figure(400 + nPage);
    subplot(nFigs_page, 1, subp);
    %%
    csCombPL_one(arrQM, prn, obs_PR, obs_CP);
    axis([tS tE -yMax yMax]);
    titleS = strcat('PRN : ', num2str(prn));
    title(titleS);
end

