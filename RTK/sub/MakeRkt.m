function R = MakeRkt(SatsInfo_Rv)

%% PR/CP의 관측치잡음(measurement noise) 초기 설정 - 고도각 고려 안함
cp1_noise = 0.01; %: Phase on L1
%% 고도각이 낮을 경우 관측치잡음을 증가 시킴
el_cut = 30;
for i=1:length(SatsInfo_Rv(:,1))
    el = SatsInfo_Rv(i,7);
    if el > el_cut
        cp1_noise_el = cp1_noise;
    else
        sin_el_2 = sind(el)*sind(el);
        cp1_noise_el = cp1_noise/sin_el_2;
%         cp1_noise_el = cosd(el)^10;
    end
    R(i,i) = cp1_noise_el;
end