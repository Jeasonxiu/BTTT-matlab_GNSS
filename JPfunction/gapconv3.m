function [target] = gapconv3(VRS_file, gap)
%
%   input VRS_file : VRS text file 

%   input gap : distance between reference receiver and target receiver
%                   ex) 43 (centi-meter)
%
%   output target : N X 4 matric = [gs, x, y, z]
%   output file : VRS_file_gap.txt

% clear all; close all;
% gap = 50;
% load('SOND_170316_adm.txt');
% VRS_file = 'SOND_170316_adm.txt';
% VRS_file ='PTCO4_hyunu_170216_adm.txt'
load(VRS_file);
VRS = eval(VRS_file(1:end-4));              
VRS_original = eval(VRS_file(1:end-4));
VRS(:,2:4) = VRS(:,5:7);
fid_out = fopen(strcat(VRS_file(1:end-4),'_',num2str(gap),'cm.txt'),'w');
gap = gap/100;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 이동시점에서 계산된 GAP만큼 shift된 좌표를 이동시점 이전까지 적용하는 과정
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
% gs(1,1) = gs(2,1) -1;
gs = VRS(:,1);
target(:,1) = gs;
for i=1:length(VRS_original(:,1))
    epoch = VRS_original(i,1);
    temp = VRS_original(find(VRS_original(:,1) == epoch),:);
    VRS_ori_line = temp(1,:);
    temp = target(find(target(:,1) == epoch),2:4);
    VRS_ori_line(5:7) = temp(1,:);
    fprintf(fid_out,'%8.3f %3.8f %3.8f %4.8f %10.8f %10.8f %10.8f %1d \n', VRS_ori_line(1),...
        VRS_ori_line(2),VRS_ori_line(3),VRS_ori_line(4),VRS_ori_line(5),VRS_ori_line(6),VRS_ori_line(7),...
        VRS_ori_line(8));
end