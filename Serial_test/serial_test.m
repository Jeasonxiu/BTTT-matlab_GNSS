% a=serial('com7','BaudRate',115200,'DataBits',8,'Parity','none','StopBits',1)
clc;
clear all;
close all;

s2=serial('com3','BaudRate',9600,'DataBits',8,'Parity','none','StopBits',1)
fopen(s2);
fid=fopen('test_test2.txt','w');
n=0;
time =1;
datalen = 100;
data = zeros(1,datalen);
% for i = 1:datalen;
while (1)
    time = time + 1;
%     data1 = query(s2, '*IDN?');
    data1 = fscanf(s2);
    disp(data1)
%     a = str2num(fscanf(s2));
%     
%     data(1:end-1) = data(2:end);
%     data(end) = a;
%     idx = time-datalen+1:1:time;
%     
% %     plot(idx, data, '*');
%     plot(idx, data, 'linewidth', 2);
%     axis([min(idx) max(idx) 0 2000]);
%     drawnow; % Plot을 계속 업데이트해줍니다. 마치 리얼타임처럼...
end

% for n = 0:5
%     n=n+1
% aaaa=fread(a);
% 
% 
% aaa=[data_Combination(aaaa)];
% 
% 
% fprintf(fid,'%s',aaa)
% end
% delete(a)