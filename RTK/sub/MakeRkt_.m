function R = MakeRkt_(SatsInfo_Rv)

for i=1:length(SatsInfo_Rv(:,1))
%     R(i,i) = 1/sind(SatsInfo_Rv(i,7));
    R(i,i) = cosd(SatsInfo_Rv(i,7))^10;
end