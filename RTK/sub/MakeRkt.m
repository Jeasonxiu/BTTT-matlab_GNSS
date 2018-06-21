function R = MakeRkt(SatsInfo_Rv)

%% PR/CP�� ����ġ����(measurement noise) �ʱ� ���� - ���� ��� ����
cp1_noise = 0.01; %: Phase on L1
%% ������ ���� ��� ����ġ������ ���� ��Ŵ
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