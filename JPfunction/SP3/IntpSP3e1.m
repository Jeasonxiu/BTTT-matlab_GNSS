function [out] = IntpSP3e1(arrSP3, prn, tt) 
%
%function [out] = IntpSP3e1(arrSP3, prn, tt)
%
% DO: Extract time-tag and interpolated coordinates X, Y, and Z for 
%     just one epoch (referred to as 'e1')
%
% <input>   arrSP3: SP3 array
%           prn: PRN number
%           tt: Time-Tag to interpolate
%
% <out>     out: Interpolated X, Y, Z
%               c1(X), c2(Y), c3(Z)
%
% Copyright: Kwan-Dong Park, December 2013 @LDEO
%
% ����: IntpSP3e1�� �� epoch�� ���ؼ�, IntpSP3�� �ܹ��� ��� epoch�� ���� ����
% -Modifications-
% 1/20/2014: tt�� �־��� �ð�����(0��) ������ �� ������ ���ܼ�, ���ǹ��� �߰��� tint<TT(1)
%           - �׷���, �� ����� ��ȣ���޽ð� ���� ���� ������ �ð����� ����ϴ� ���� ���� ����
% 11/30/2014: ���out�� row ���ͷ� ����
% 8/8/2015: ��� out���� �ð� ���� - Ÿ �Լ��� ���� �����˵��� ���! ���� ���� _old

%% �Էº��� �׽�Ʈ
% load 'SP3'; arrSP3 = sp3;
% prn = 1;
% pos = 1;
% tt =  55;
%% �Էµ� arrSP3�� TIME-TAG ���� �� ���� ����, ���̴� 96�̾�� ��
TT = unique(arrSP3(:,1));
nTT = length(TT);
if nTT ~= 96, error('The number of epochs is not equal to 96'); end
%% �ش� PRN�� ��ǥ ����: posArr�� 3�̸� X, 4�̸� Y, 5�̸� Z, 6�̸� dt(�ð����)
id1 = find(arrSP3(:, 2) == prn);
if isempty(id1)
    out = [];
    return
end
    
tPRN = arrSP3(id1, 1);
xPRN = arrSP3(id1, 3); %: 3=X, 4=Y, 5=Z, 6=dt
yPRN = arrSP3(id1, 4);
zPRN = arrSP3(id1, 5); 
%% Lagrange Interpolation�� ���� ������ 9�� ���õǾ� ����
nOrder = 9; %: ����μ� �ݵ�� Ȧ���� ����ؾ� ��. nRange�� �Ҽ��� ���� ��
nPoints = nOrder - 1;
nRange = (nOrder - 1) / 2;
%%
out = zeros(1, 3);  %: ������ ũ�� ���� ����
tint = tt;          %: �����Ϸ��� �ð�
iTT = PickTTrange(TT, tint);

if tint < TT(1) || iTT - nRange <= 0 %: �ð迭 ���� �ð��� �Էµ� ��쿡 ���
    iStart = 1;
    iEnd = nPoints;
elseif iTT + nRange > nTT
    iStart = nTT - nPoints + 1;
    iEnd = nTT;
else
    iStart = iTT - nRange + 1;
    iEnd = iTT + nRange;
end
%%
tin = tPRN(iStart:iEnd);
xin = xPRN(iStart:iEnd);
yin = yPRN(iStart:iEnd);
zin = zPRN(iStart:iEnd);
xout = lagrange(tin, xin, tint);
yout = lagrange(tin, yin, tint);
zout = lagrange(tin, zin, tint);
out(1) = xout;
out(2) = yout;
out(3) = zout;
