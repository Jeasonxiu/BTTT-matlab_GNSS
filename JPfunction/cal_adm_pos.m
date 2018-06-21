function [adm_pos]=cal_adm_pos(pos,time)
% cal_adm_pos
% 연간(Annual) 변위량(Displacement)을 반영(Modified)한 좌표(POS) 산출
% SUWN 상시관측소 속도 기반
% Dusik Kim, 2015.06.29
% pos=[x y z];
% time=[doy year];
%-Input File Test-----------------------------
% pos=[-3026443.3260  4067351.6671  3857212.8707];
% time=[329, 2014];
%---------------------------------------------

t0=yrscale(001, 2000);
t=yrscale(time(1), time(2));

vx=-0.270777095438770E-01;
vy=-0.101986772893401E-01;
vz=-0.946613643617538E-02;

% 속도를 반영하여 산출한 XYZ 좌표
adm_pos=[pos(1,1)+vx*(t-t0) pos(1,2)+vy*(t-t0) pos(1,3)+vz*(t-t0)];
