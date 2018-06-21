function RTKlite_postError(GNGGA, result);
tHour = [1:1:length(result(:,1))]';

titlestr = ({['Horizontal RMSE = ', num2str(rms(result(:,4)))];...
    ['3D RMSE = ', num2str(rms(result(:,5)))]});

figure();
hold on;

subplot(4,4,[1,2,5,6])
% xlim([-1,1]);
% ylim([-1,1]);
% if ~isempty(
plot(result(find(result(:,6) == 4),2), result(find(result(:,6) == 4),1),'bo'); hold on;
plot(result(find(result(:,6) == 5),2), result(find(result(:,6) == 5),1),'r*')
plot(result(find(result(:,6) == 1),2), result(find(result(:,6) == 1),1),'k^'); 
% lineCross(0,0,'r',0.1)
% legend('Not Fix','Float','Fix')
title(titlestr)
axis([-1 1 -1 1]);
axis([-max(max(abs(result(:,1:2)))) max(max(abs(result(:,1:2))))...
    -max(max(abs(result(:,1:2)))) max(max(abs(result(:,1:2))))]);
axis square
grid on; 
% xlabel({'\Delta E (meters)',...
%     ['dNE RMSE = ', num2str(decimal(rms(result(:,4)))), '   std =', num2str(decimal(std(result(:,4))))],...
%     [' 3D RMSE = ', num2str(decimal(rms(result(:,5)))), '   std =', num2str(decimal(std(result(:,5))))]}); 
% xlabel({'\Delta E (meters)',...
%     ['dNE RMSE = ', num2str(decimal(rms(result(:,4)))), '   std =', num2str(decimal(std(result(:,4)))),...
%     ' 3D RMSE = ', num2str(decimal(rms(result(:,5)))), '   std =', num2str(decimal(std(result(:,5))))]}); 
ylabel('\Delta N (meters)')
legend('Fix', 'Float')

subplot(4,4,[3,4])
plot(result(find(result(:,6) == 4),10), result(find(result(:,6) == 4),1), '.r:',...
    result(find(result(:,6) == 4),10), result(find(result(:,6) == 4),2), '.b:'); xlim([1 length(tHour)]);
grid on; hold on;
plot(result(find(result(:,6) == 5),10), result(find(result(:,6) == 5),1), '*r',...
    result(find(result(:,6) == 5),10), result(find(result(:,6) == 5),2), '*b'); xlim([1 length(tHour)]);
plot(result(find(result(:,6) == 1),10), result(find(result(:,6) == 1),1), '^r',...
    result(find(result(:,6) == 1),10), result(find(result(:,6) == 1),2), '^b'); xlim([1 length(tHour)]);
ylim([-0.025 0.025]);
legend('\Delta N(Fix)', '\Delta E(Fix)','\Delta N(Float)', '\Delta E(Float)')
ylabel('\Delta N,E (meters)');

subplot(4,4,[7,8])
plot(result(find(result(:,6) == 4),10), result(find(result(:,6) == 4),3), 'b.:'); xlim([1 length(tHour)]); 
grid on; hold on;
plot(result(find(result(:,6) == 5),10), result(find(result(:,6) == 5),3), 'r*'); xlim([1 length(tHour)]); 
plot(result(find(result(:,6) == 1),10), result(find(result(:,6) == 1),3), 'k^'); xlim([1 length(tHour)]); 
ylim([-0.1 0.1]);
ylabel('\Delta V (meters)')
legend('\Delta V(Fix)','\Delta V(Float)')

subplot(4,4,[9,10])
plot(result(find(result(:,6) == 4),4), '.b:'); xlim([1 length(tHour)]); grid on; hold on;
plot(result(find(result(:,6) == 5),10), result(find(result(:,6) == 5),4), 'r*');
plot(result(find(result(:,6) == 1),10), result(find(result(:,6) == 1),4), 'k^:');
ylim([0 0.05]);
ylabel('Error (meters)')
legend('\Delta NE')

subplot(4,4,[11,12])
plot(tHour(:,1), result(:,6),'.b:'); 
xlim([1 length(tHour)]); grid on;
ylim([1 7]);
ylabel('Fix Quality');

subplot(4,4,[13,14])
plot(result(find(result(:,6) == 4),5), '.b:'); xlim([1 length(tHour)]); grid on; hold on;
plot(result(find(result(:,6) == 5),10), result(find(result(:,6) == 5),5), 'r*')
plot(result(find(result(:,6) == 1),10), result(find(result(:,6) == 1),5), 'k^'); 
ylim([0 0.1]);
ylabel('Error (meters)')
legend('\Delta 3D')

subplot(4,4,[15,16])
plot(result(find(result(:,6) == 4),10), result(find(result(:,6) == 4),8), 'b.:'); xlim([1 length(tHour)]); grid on; hold on;
plot(result(:,10), result(:,7), 'r.:'); xlim([1 length(tHour)]); 
plot(result(find(result(:,6) == 5),10), result(find(result(:,6) == 5),8), 'r*'); 
plot(result(find(result(:,6) == 1),10), result(find(result(:,6) == 1),8), 'k^'); 
legend('Num of all Sats', 'Num of GPS')
ylabel('Number of Satellites')

figure()
subplot(3,1,[1])
plot(result(find(result(:,6) == 4),4), '.b:'); xlim([1 length(tHour)]); grid on; hold on;
plot(result(find(result(:,6) == 5),10), result(find(result(:,6) == 5),4), 'r*');
plot(result(find(result(:,6) == 1),10), result(find(result(:,6) == 1),4), 'k^:');
% ylim([0 0.2]);
ylabel('Error (meters)')
subplot(3,1,[2])
plot(result(find(result(:,6) == 4),5), '.b:'); xlim([1 length(tHour)]); grid on; hold on;
plot(result(find(result(:,6) == 5),10), result(find(result(:,6) == 5),5), 'r*')
plot(result(find(result(:,6) == 1),10), result(find(result(:,6) == 1),5), 'k^'); 
% ylim([0 0.2]);
ylabel('Error (meters)')
legend('\Delta 3D')
legend('\Delta NE')
subplot(3,1,[3])
plot(GNGGA(:,14), '.b:'); xlim([1 length(tHour)]); grid on; hold on;


% 
% figure();
% hold on;
% 
% subplot(5,9,[1,2,3,4,10,11,12,13,...
%     19,20,21,22,28,29,30,31])
% xlim([-0.08,0.08]);
% ylim([-0.08,0.08]);
% plot(result(find(result(:,6) == 4),2), result(find(result(:,6) == 4),1),'b.'); hold on;
% lineCross(0,0,'r',0.1)
% % legend('Not Fix','Float','Fix')
% axis([-0.008 0.008 -0.008 0.008]);
% axis square
% grid on; 
% % xlabel({'\Delta E (meters)',...
% %     ['dNE RMSE = ', num2str(decimal(rms(result(:,4)))), '   std =', num2str(decimal(std(result(:,4))))],...
% %     [' 3D RMSE = ', num2str(decimal(rms(result(:,5)))), '   std =', num2str(decimal(std(result(:,5))))]}); 
% % xlabel({'\Delta E (meters)',...
% %     ['dNE RMSE = ', num2str(decimal(rms(result(:,4)))), '   std =', num2str(decimal(std(result(:,4)))),...
% %     ' 3D RMSE = ', num2str(decimal(rms(result(:,5)))), '   std =', num2str(decimal(std(result(:,5))))]}); 
% ylabel('\Delta N (meters)')
% xlabel('\Delta E (meters)')
% 
% 
% subplot(5,9,[6,7,8,9])
% plot(result(find(result(:,6) == 4),10), result(find(result(:,6) == 4),1), '.b'); xlim([1 length(tHour)]);
% grid on; hold on;
% ylim([-0.008 0.008]);
% ylabel('\Delta N(meters)')
% 
% 
% subplot(5,9,[15,16,17,18])
% plot(result(find(result(:,6) == 4),10), result(find(result(:,6) == 4),2), 'b.'); xlim([1 length(tHour)]); 
% grid on; hold on;
% ylim([-0.008 0.008]);
% ylabel('\Delta E(meters)')
% 
% 
% subplot(5,9,[33,34,35,36])
% plot(result(find(result(:,6) == 4),4), '.b'); xlim([1 length(tHour)]); grid on; hold on;
% ylim([0 0.008]);
% ylabel('H Error(meters)')
% 
% 
% subplot(5,9,[24,25,26,27])
% plot(result(find(result(:,6) == 4),10), result(find(result(:,6) == 4),3), 'b.'); xlim([1 length(tHour)]); 
% grid on; hold on;
% ylim([-0.02 0.02]);
% ylabel('\Delta V(meters)')
% 
% 
% subplot(5,9,[42,43,44,45])
% plot(result(find(result(:,6) == 4),5), 'b.'); xlim([1 length(tHour)]); grid on; hold on;
% ylim([0 0.02]);
% ylabel('3D Error(meters)')
% 
% 
% 
% figure();
% hold on;
% xlim([-0.08,0.08]);
% ylim([-0.08,0.08]);
% plot(result(find(result(:,6) == 4),2), result(find(result(:,6) == 4),1),'b.'); hold on;
% lineCross(0,0,'r',0.1)
% % lineCross(0,0,'r',0.1)
% % legend('Not Fix','Float','Fix')
% axis([-0.008 0.008 -0.008 0.008]);
% axis square
% grid on; 
% % xlabel({'\Delta E (meters)',...
% %     ['dNE RMSE = ', num2str(decimal(rms(result(:,4)))), '   std =', num2str(decimal(std(result(:,4))))],...
% %     [' 3D RMSE = ', num2str(decimal(rms(result(:,5)))), '   std =', num2str(decimal(std(result(:,5))))]}); 
% % xlabel({'\Delta E (meters)',...
% %     ['dNE RMSE = ', num2str(decimal(rms(result(:,4)))), '   std =', num2str(decimal(std(result(:,4)))),...
% %     ' 3D RMSE = ', num2str(decimal(rms(result(:,5)))), '   std =', num2str(decimal(std(result(:,5))))]}); 
% ylabel('\Delta N (meters)')
% xlabel('\Delta E (meters)')