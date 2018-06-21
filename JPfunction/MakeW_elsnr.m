function [W] = MakeW_elsnr(el,S1)
%
%function [W] = MakeW_elpr(el,S1)
%
% DO: Assign an weight for the given elevation angle [pseudo-range]
%
% <input>   el: elevation angle
%
% <output>  W: Weight 
%
% Copyright: Kwan-Dong Park, Jipyong Space, October 18, 2014
%

%% 코드-의사거리 가중치 기본은 1미터로 설정
pr_weight = 1.00;
%% 일정 고도각(el_cut) 이상은 동일한 가중치 부여
% el_cut = 30;
% if el > el_cut
%     W = 1  ;
% else
%     W = pr_weight*cscd(el);
% end
%% Divide by Zero 방지 - 1도 이하일때 1도 값으로 설정
if el < 1
    W_el = pr_weight*cscd(1);
else
    W_el = pr_weight*cscd(el);
end
% if el > 50 
%     W_el = 1;
% else
%     W_el = pr_weight*cscd(el);
% end
%% 역치로 리턴하기
W_el = (1/W_el)^2;
% SNR = S1;
SNR = 2.225 * (S1) - 32;
if SNR > 80
SNR = 90;
end
% W_snr = sind(SNR);
W_snr = sind(SNR)^2;
% W = W_el * W_snr;
W = W_snr;

%% 10/18/2014 김혜인 고도각 가중치
% W = sind(el)^2;