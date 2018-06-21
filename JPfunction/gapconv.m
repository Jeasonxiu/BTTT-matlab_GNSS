function [target] = gapconv(VRS, gap)
%
%   input VRS : VRS [N X 4] matrix -> gs, x ,y ,z

%   input gap : distance between reference receiver and target receiver
%                   ex) 0.43 (meter)
%
%   output target : N X 4 matric = [gs, x, y, z]

% clear all; close all;
% gap = 0.50;
% load('VRS_170216_set1.mat');
% VRS = vrs_joon;
% VRS(:,2:4) = VRS(:,5:7);

for i =2:length(VRS(:,1))
    
    A = VRS(i,2:4);
    B = VRS(i-1,2:4);
    dis(i,:) = norm(A-B); Dis(i,1) = VRS(i,1); Dis(i,2) = dis(i,1);
    target(i,1) = VRS(i,1);
    target(i,2:4) = [0,0,0]; C = [0, 0, 0];
    gs(i,1) = VRS(i,1);

    if dis(i) > 0.1                 % VRS 움직임 시작
        m = (dis(i) - gap) / dis(i); M(i,1) = m;
        n = gap / dis(i); M(i,2) = n;
        x = (m * A(1) + n * B(1)) ;
        y = (m * A(2) + n * B(2)) ;
        z = (m * A(3) + n * B(3)) ;
        target(i,1) = VRS(i,1);
        target(i,2:4) = [x,y,z]; C = [x, y, z];

    elseif dis(i) < gap && target(i-1,2) ~= 0       % target의 이전좌표가 존재하나 VRS 이동거리가 gap보다 작을때 VRS 차이값 만큼 target 좌표 생성
        AA(i-1,1) = 1;
        C = target(i-1,2) + A - B;
        target(i,2:4) = target(i-1,2:4) + A - B;

    end

end
% gs = round(gs);
i = 1;
udVRS = flipud(VRS);
udtarget = flipud(target);
while target(i,2) == 0
    i = i + 1;
    start = target(i,1);
end
for j = find(udVRS(:,1) == start) : length(VRS)-1
    invT = udVRS(j,1);
    A = udVRS(j,2:4);
    B = udVRS(j+1,2:4);
    C = A - B;
    udtarget(j+1,2:4) = udtarget(j,2:4) - C;
    udtarget(j+1,1) = invT;
end
target = flipud(udtarget);
gs(1,1) = gs(2,1) -1;
target(:,1) = gs;