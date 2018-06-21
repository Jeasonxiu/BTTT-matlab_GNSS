function N_before = makeNbefore(N, prn_RS, SatsInfo_Rv);
k = 0;
for i = 1:32
    prn = i;
    if (SatsInfo_Rv(i,1) == prn) & (SatsInfo_Rv(i,1) ~= prn_RS)
        k = k + 1;
        N_before(i,1) = N(k);
    else
        N_before(i,1) = 0;
    end
end