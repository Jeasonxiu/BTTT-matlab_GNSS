function [hour] = gs2h24(gs)

second = mod(gs,86400);
hour = second/3600;




end


