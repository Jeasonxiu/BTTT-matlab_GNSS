function [xyz,vxyz,xyzLS] = RK4_my(x, deltat,fid_new)
% 
% function [xyz,vxyz] = RK4(x, deltat)
%                      Runge-Kutta 4차 적분
% <input>:
%   x = initial value          ex) [1prn 2gs 3x 4y 5z 6xdot 7ydot 8zdot 9xLS 10yLS 11zLS ...]
%   deltat= integration term   ex) 10
%   fid_new 
% 
% <output>:
%   xyz = RK-4차를 적용하여 산출한 GLONASS의 좌표
%   vxyz = xyz 위치에서의 GLONASS의 속도; 상대성 효과 보정에 사용하기 위함
% 
%  Copyright: JongSeok Kim, July 1, 2016

%%
x0 = x; % 초기값 저장
acc = x(9:11);
x=x(3:8);
dt = deltat; 

k1= EOM_glo_my(x,acc); % GLONASS 운동방정식, 입출력만 변경 안에 내용은 동일
w1=x + (k1.*dt)/2.0;
k2 = EOM_glo_my(w1,acc);
w2=x + (k2.*dt)/2.0;
k3 = EOM_glo_my(w2,acc);
w3=x + (k3.*dt);
k4 = EOM_glo_my(w3,acc);
x = x+ (k1 + (2.0.*k2) + (2.0.*k3) + k4) *dt / 6.0;
xyz=x(1:3);
vxyz=x(4:6);
xyzLS = acc;
% fprintf(fid_my,'gs : %8.3f %8.3f %8.3f\n',x0(1), x0(2),dt);
% fprintf(fid_my,'x0 : %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',x0(3),x0(4),x0(5),x0(6),x0(7),x0(8),x0(9),x0(10),x0(11) );
% fprintf(fid_my,'k1 : %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n',k1(1),k1(2),k1(3),k1(4),k1(5),k1(6));
% fprintf(fid_my,'w1 : %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n',w1(1),w1(2),w1(3),w1(4),w1(5),w1(6));
% fprintf(fid_my,'k2 : %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n',k2(1),k2(2),k2(3),k2(4),k2(5),k2(6));
% fprintf(fid_my,'w2 : %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',w2(1),w2(2),w2(3),w2(4),w2(5),w2(6) );
% fprintf(fid_my,'k3 : %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n',k3(1),k3(2),k3(3),k3(4),k3(5),k3(6));
% fprintf(fid_my,'w3 : %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n',w3(1),w3(2),w3(3),w3(4),w3(5),w3(6) );
% fprintf(fid_my,'k4 : %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n',k4(1),k4(2),k4(3),k4(4),k4(5),k4(6));
% fprintf(fid_my,'x : %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n',x(1),x(2),x(3),x(4),x(5),x(6));
% fprintf(fid_my,'------------------------------------------------------------\n');
