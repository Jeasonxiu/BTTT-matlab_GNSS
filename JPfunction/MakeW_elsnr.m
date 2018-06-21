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

%% �ڵ�-�ǻ�Ÿ� ����ġ �⺻�� 1���ͷ� ����
pr_weight = 1.00;
%% ���� ����(el_cut) �̻��� ������ ����ġ �ο�
% el_cut = 30;
% if el > el_cut
%     W = 1  ;
% else
%     W = pr_weight*cscd(el);
% end
%% Divide by Zero ���� - 1�� �����϶� 1�� ������ ����
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
%% ��ġ�� �����ϱ�
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

%% 10/18/2014 ������ ���� ����ġ
% W = sind(el)^2;