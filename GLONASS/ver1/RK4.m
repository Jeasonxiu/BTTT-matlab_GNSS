function [xyz,vxyz] = RK4(x, deltat)
% 
% function [xyz,vxyz] = RK4(x, deltat)
%                      Runge-Kutta 4�� ����
% <input>:
%   x = initial value          ex) [prn gs x y z xdot ydot zdot xLS yLS zLS ...]
%   deltat= integration term   ex) 10
%
% <output>:
%   xyz = RK-4���� �����Ͽ� ������ GLONASS�� ��ǥ
%   vxyz = xyz ��ġ������ GLONASS�� �ӵ�; ��뼺 ȿ�� ������ ����ϱ� ����
% 
%  Copyright: Mi-so Kim, January 14th, 2015

%% �Լ� ������� test
% eph_file ='brdc2120.14g'; [eph_glo] = ReadEPH_GLO(eph_file);
% icol = 3; x = eph_glo(icol,:);
% deltat = 10;
 
%%
x0 = x; % �ʱⰪ ����
N = 4;  % ���� ����; Runge-Kutta 4th
dt  = deltat; dt2=dt*dt; 
ddt1=[2 8; 2 8; 1 2;]; % ����� Ȯ�� case1 - ȣ������ �ڵ�
% ddt1=[2 4; 2 4; 1 2;]; % case 2 - �ܼ������̶�� case 2�� �´µ�..?
for i=1:N
    f(i,:) = EOM_glo(x(1,:)); % GLONASS �������
    if i == N
        xyz = x0(3:5) + x0(6:8)*dt + (f(1,:)+2*(f(2,:)+f(3,:))+f(4,:))*dt2/12; % ��ġ �� �ӵ� ���: 4��
        vxyz = x0(6:8) + (f(1,:)+2*(f(2,:)+f(3,:))+f(4,:))*dt/6; 
    else
        p_n(i,:) = x0(3:5) + x0(6:8)*(dt/ddt1(i,1)) + f(i,:)*(dt2/ddt1(i,2)); % ��ġ �� �ӵ� ���: 1-3��
        v_n(i,:) = x0(6:7) + f(i,1:2)*(dt/ddt1(i,1));
        x(3:5) = p_n(i,:); x(6:7) = v_n(i,:);
    end
end

% x0�� �����ϴ� ���� �³�? / �����Ÿ�� ���� ������ ��ǥ�� ����?