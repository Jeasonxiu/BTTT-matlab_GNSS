function [R] = MakeR_rtap_l1(R, indxRR, el)
% % %
% % %function [R] = MakeR_rtap_sf(R, indxRR, el)
% % %
% % % DO: Create a measurement noise matrix 
% % %
% % % <input>   R: measurement noise matrix
% % %           indxRR: row number to place the measurement noise
% % %           el: elevation angle [degrees]
% % %
% % % <output>  R: measurement noise matrix
% % %
% % % Copyright: Kwan-Dong Park @LDEO, March 31, 2014
% % % --- Modifications ---
% % % 4/2/14 cosine�� �ƴ� sine���� ��ü
% % % 1/7/2015 ���� ���ļ� ������� ���� (Jihye Won)

%% PR/CP�� ����ġ����(measurement noise) �ʱ� ���� - ���� ��� ����
cp1_noise = 0.01; %: Phase on L1
%% ������ ���� ��� ����ġ������ ���� ��Ŵ
el_cut = 30;
if el > el_cut
    cp1_noise_el = cp1_noise;
else
    sin_el_2 = sind(el)*sind(el);
    cp1_noise_el = cp1_noise/sin_el_2;
end
%% ����ġ ���� ���
R(indxRR, indxRR) = cp1_noise_el;
