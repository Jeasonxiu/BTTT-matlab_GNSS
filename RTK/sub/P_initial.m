function P = P_initial(SatsInfo, SatsInfo_before, prn_RS, P, P_before)

Sats_compare = SatsInfo - SatsInfo_before;
Sats_compare = Sats_compare(find(Sats_compare ~= prn_RS));
riselist = find(Sats_compare > 0);
setlist = find(Sats_compare < 0);
%% P 炼钦 积己
if length(riselist) > 0
    for i = 1:length(riselist)
        prn_OS = riselist(i);
        P_before(prn_OS,1) = 100;
    end
end
%% RS row 客 set 困己 昏力
if length(setlist) > 0
    P_before(setlist) = 0;
    P_before(prn_RS) = 0;
end
P_temp = P_before(find(P_before(:,1) ~= 0),:);

for i = 1:length(P_temp)
    P(i+3,i+3) = P_temp(i);
end