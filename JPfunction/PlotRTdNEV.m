function [] = PlotRTdNEV(TruePos, EstPos, kE, NoEpochs)
%
%function [] = PlotRTdNEV(TruePos, EstPos, kE, NoEpochs)
%
% DO: Plot real-time positional error 
%
% <input>   TruePos: True position in XYZ
%           EstPos: Estimated position in XYZ
%           kE: Epoch index (1-second sampling assumed)
%           NoEpochs: Number of epochs to process
%
% <output>  Real-time Plot
%
% Copyright: Kwan-Dong Park, 1/8/2015 @Jipyong Space
% Modified by Jihye Won, 1/9/2015 @Jipyong Space

%% �׷����� ǥ���� �ð� ����
div = 60; %: �д��� ǥ��
tS = 0;
tE = NoEpochs/div;
tt = kE/div;
%% ������ �ش��ϴ� ���浵 ���� ����
gd = xyz2gd(TruePos); TrueLat = gd(1); TrueLon = gd(2); 
%% ������ XYZ�� ���� XYZ�� ���̰� ����
dXYZ = EstPos - TruePos;
%% dXYZ�� dNEV�� ��ȯ
dNEV = xyz2topo(dXYZ, TrueLat, TrueLon);
dN = dNEV(1); dE = dNEV(2); dV = dNEV(3);
hst = num2str(sqrt(dE^2+dN^2),'%4.2f');
hsl = strcat('H error:',hst,'(m)');
dst = num2str(sqrt(dV^2),'%4.2f');
dsl = strcat('V error:',dst,'(m)');

%% �׷��� 111���� �� ���� ���� ǥ��
figure(111);
set(gcf,'position',[775    54   580   631]);
%% ���� ����
subplot(4,1,1:3);
plot([0 0],[-5 5],'k-','linewidth',1);
plot([-5 5],[0 0],'k-','linewidth',1);
plot(dE, dN, 'ob','markerfacecolor',[1-(1/NoEpochs) 0 1/NoEpochs]); 
% if norm(dNEV(1:2)) < 0.3; 
%     axis([-0.5 0.5 -0.5 0.5]);
if norm(dNEV(1:2)) < 1.5; 
    axis([-1.5 1.5 -1.5 1.5]);
% elseif norm(dNEV(1:2)) < 0.5;    
%     axis([-1 1 -1 1]);
else
    axis([-5 5 -5 5]); 
end
grid on; hold on
% xlabel('\Delta E(m)');
ylabel('\Delta N (m)');
title(hsl);
%% ���� ����
subplot(4,1,4);
plot(tt, dV, 'ob','markerfacecolor',[1-(1/NoEpochs) 0 1/NoEpochs]); 
% legend(dsl);
% if abs(dNEV(3)) < 0.5;
%     axis([tS tE -0.75 0.75]);
if abs(dNEV(3)) < 1.5;
    axis([tS tE -2 2]);
elseif abs(dNEV(3)) < 3;
    axis([tS tE -5 5]);
else
    axis([tS tE -10 10]); 
end
grid on; hold on
xlabel('Time Elapsed (Minutes)');
ylabel('\Delta V (m)');
title({'\Delta E (m)',dsl})