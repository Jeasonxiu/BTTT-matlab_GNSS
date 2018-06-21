function P_before = makePbefore(P, prn_RS, SatsInfo_Rv);
k = 1;
for i = 1:32
    prn = i;
    if (SatsInfo_Rv(i,1) == prn) & (SatsInfo_Rv(i,1) ~= prn_RS)
        P_before(i,1) = P(k+3,k+3);
        k = k + 1;
    else
        P_before(i,1) = 0;
    end
end

