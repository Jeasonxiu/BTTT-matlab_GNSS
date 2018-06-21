function R = MakeRrtk_Pr(SatsList_OS, SatsInfo_Rv,QMRv_Pr_1e)
prnlist = intersect(SatsList_OS, SatsInfo_Rv(:,1));

for i=1:length(prnlist(:,1))
    prn_OS = prnlist(i);
    Pr = QMRv_Pr_1e(find(QMRv_Pr_1e(:,2) == prn_OS),4);
    R(i,i) = Pr;
end