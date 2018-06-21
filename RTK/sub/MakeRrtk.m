function R = MakeRrtk(SatsList_OS, SatsInfo_Rv)
prnlist = intersect(SatsList_OS, SatsInfo_Rv(:,1));
cp1_noise = 0.01; %: Phase on L1
for i=1:length(prnlist(:,1))
    prn_OS = prnlist(i);
    el = SatsInfo_Rv(find(SatsInfo_Rv(:,1) == prn_OS),7);
    if el > 30
        R(i,i) = cp1_noise;
    else
        R(i,i) = cp1_noise/sind(SatsInfo_Rv(find(SatsInfo_Rv(:,1) == prn_OS),7));
% R(i,i) = cosd(SatsInfo_Rv(find(SatsInfo_Rv(:,1) == prn_OS),7))^10;
    end
end