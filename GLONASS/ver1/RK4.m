function [xyz,vxyz] = RK4(x, deltat)
% 
% function [xyz,vxyz] = RK4(x, deltat)
%                      Runge-Kutta 4차 적분
% <input>:
%   x = initial value          ex) [prn gs x y z xdot ydot zdot xLS yLS zLS ...]
%   deltat= integration term   ex) 10
%
% <output>:
%   xyz = RK-4차를 적용하여 산출한 GLONASS의 좌표
%   vxyz = xyz 위치에서의 GLONASS의 속도; 상대성 효과 보정에 사용하기 위함
% 
%  Copyright: Mi-so Kim, January 14th, 2015

%% 함수 만들기전 test
% eph_file ='brdc2120.14g'; [eph_glo] = ReadEPH_GLO(eph_file);
% icol = 3; x = eph_glo(icol,:);
% deltat = 10;
 
%%
x0 = x; % 초기값 저장
N = 4;  % 적분 차수; Runge-Kutta 4th
dt  = deltat; dt2=dt*dt; 
ddt1=[2 8; 2 8; 1 2;]; % 상수값 확인 case1 - 호석오빠 코드
% ddt1=[2 4; 2 4; 1 2;]; % case 2 - 단순적분이라면 case 2가 맞는듯..?
for i=1:N
    f(i,:) = EOM_glo(x(1,:)); % GLONASS 운동방정식
    if i == N
        xyz = x0(3:5) + x0(6:8)*dt + (f(1,:)+2*(f(2,:)+f(3,:))+f(4,:))*dt2/12; % 위치 및 속도 계산: 4차
        vxyz = x0(6:8) + (f(1,:)+2*(f(2,:)+f(3,:))+f(4,:))*dt/6; 
    else
        p_n(i,:) = x0(3:5) + x0(6:8)*(dt/ddt1(i,1)) + f(i,:)*(dt2/ddt1(i,2)); % 위치 및 속도 계산: 1-3차
        v_n(i,:) = x0(6:7) + f(i,1:2)*(dt/ddt1(i,1));
        x(3:5) = p_n(i,:); x(6:7) = v_n(i,:);
    end
end

% x0에 갱신하는 것이 맞나? / 룬게쿠타로 계산된 이전의 좌표에 갱신?