function [W] = MakeW_elpr(el)
%
%function [W] = MakeW_elpr(el)
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
    W = pr_weight*cscd(1);
else
    W = pr_weight*cscd(el);
end
%% ��ġ�� �����ϱ�
% W = (1/W);
W = (1/W)^2;
%% 10/18/2014 ������ ���� ����ġ
% W = sind(el)^2;