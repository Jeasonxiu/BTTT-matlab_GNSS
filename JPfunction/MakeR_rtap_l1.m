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
% % % 4/2/14 cosine이 아닌 sine으로 대체
% % % 1/7/2015 단일 주파수 기반으로 수정 (Jihye Won)

%% PR/CP의 관측치잡음(measurement noise) 초기 설정 - 고도각 고려 안함
cp1_noise = 0.01; %: Phase on L1
%% 고도각이 낮을 경우 관측치잡음을 증가 시킴
el_cut = 30;
if el > el_cut
    cp1_noise_el = cp1_noise;
else
    sin_el_2 = sind(el)*sind(el);
    cp1_noise_el = cp1_noise/sin_el_2;
end
%% 관측치 잡음 행렬
R(indxRR, indxRR) = cp1_noise_el;
