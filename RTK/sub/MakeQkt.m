function Q = MakeQkt(NoSats)
Q = eye(3+NoSats);
% Q(1:3,1:3) = Q(1:3,1:3) * 1e-5;                 % Q xyz �ʱⰪ ����
% Q(4:3+NoSats,4:3+NoSats) = Q(4:3+NoSats,4:3+NoSats) * 10e-8;    % Q ���� �� N �ʱⰪ ����
Q(1:3,1:3) = Q(1:3,1:3)*10;
Q(4:3+NoSats,4:3+NoSats) = Q(4:3+NoSats,4:3+NoSats) * 0.001;    % Q ���� �� N �ʱⰪ ����

% Q(1:3,1:3) = Q(1:3,1:3)*0.03^2;
% Q(4:3+NoSats,4:3+NoSats) = Q(4:3+NoSats,4:3+NoSats) * 0.0003^2;    % Q ���� �� N �ʱⰪ ����







% %% 180515 S8 parameter
% Q(1:3,1:3) = Q(1:3,1:3)*0.19^2;
% Q(4:3+NoSats,4:3+NoSats) = Q(4:3+NoSats,4:3+NoSats) * 6;    % Q ���� �� N �ʱⰪ ����


