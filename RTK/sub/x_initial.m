function x = x_initial(SatsInfo, SatsInfo_before, prn_RS, x, x_before)

Sats_compare = SatsInfo - SatsInfo_before;
Sats_compare = Sats_compare(find(Sats_compare ~= prn_RS));
riselist = find(Sats_compare > 0);
setlist = find(Sats_compare < 0);
%% x 炼钦 积己
if length(riselist) > 0
    for i = 1:length(riselist)
        prn_OS = riselist(i);
        x_before(prn_OS,1) = 500;
    end
end
%% RS row 客 set 困己 昏力
if length(setlist) > 0
    x_before(setlist) = 0;
    x_before(prn_RS) = 0;
end
x_temp = x_before(find(x_before(:,1) ~= 0),:);
x_new = x(1:3);
for i = 1:length(x_temp)
    x_new(i+3) = x_temp(i);
end
x = x_new;