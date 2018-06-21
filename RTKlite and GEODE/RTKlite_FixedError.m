function RTKlite_FixedError(Fixed, result_fixed);
tHour = mod(Fixed(:,1), 86400);
tHour = tHour/3600;
nextdayidx = find(Fixed(:,2)==23 & Fixed(:,3) == 59 & Fixed(:,4) == 59);
if ~isempty(nextdayidx)
    tHour(nextdayidx+1:length(tHour),1) = tHour(nextdayidx+1:length(tHour),1) +24;
end
titlestr = ({['Horizontal RMSE = ', num2str(rms(result_fixed(:,4)))];...
    ['3D RMSE = ', num2str(rms(result_fixed(:,5)))]});

figure();
hold on;
subplot(4,4,[1,2,5,6])
plot(result_fixed(:,2), result_fixed(:,1),'bo'); hold on;
title(titlestr)
% axis([-0.5 0.5 -0.5 0.5]);
axis([-max(max(abs(result_fixed(:,1:2)))) max(max(abs(result_fixed(:,1:2))))...
    -max(max(abs(result_fixed(:,1:2)))) max(max(abs(result_fixed(:,1:2))))]);
axis square
grid on; 
ylabel('\Delta H (meters)')

subplot(4,4,[3,4])
plot(tHour, result_fixed(:,1), '.b:');
xlim([min(tHour) max(tHour)]);
grid on; hold on;
ylim([-max(abs(result_fixed(:,1))) max(abs(result_fixed(:,1)))]);
ylabel('\Delta N (meters)');
xlabel('Hours');

subplot(4,4,[7,8])
plot(tHour, result_fixed(:,2), '.b:');
xlim([min(tHour) max(tHour)]);
grid on; hold on;
ylim([-max(abs(result_fixed(:,2))) max(abs(result_fixed(:,2)))]);
ylabel('\Delta E (meters)');
xlabel('Hours');

subplot(4,4,[9,10])
plot(tHour, result_fixed(:,4), '.b:');
grid on; hold on;
ylim([0 max(abs(result_fixed(:,4)))]);
xlim([min(tHour) max(tHour)]);
ylabel('Horizintal Error (meters)')
xlabel('Hours');

subplot(4,4,[11,12])
plot(tHour, result_fixed(:,3), '.b:');
grid on; hold on;
ylim([-max(abs(result_fixed(:,3))) max(abs(result_fixed(:,3)))]);
xlim([min(tHour) max(tHour)]);
ylabel('\Delta V (meters)')
xlabel('Hours');

subplot(4,4,[13,14])
plot(tHour, result_fixed(:,5), '.b:');
grid on; hold on;
ylim([0 max(abs(result_fixed(:,5)))]);
xlim([min(tHour) max(tHour)]);
ylabel('3D Error (meters)')
xlabel('Hours');

subplot(4,4,[15,16])
plot(tHour, result_fixed(:,8), '.b:');
grid on; hold on;
% ylim([0 0.05]);
xlim([min(tHour) max(tHour)]); 
ylabel('Number of Satellites')

figure()
subplot(3,1,[1])
plot(tHour, result_fixed(:,4), '.b:');
grid on; hold on;
ylim([0 max(abs(result_fixed(:,4)))]);
xlim([min(tHour) max(tHour)]);
ylabel('Horizontal Error (meters)')
xlabel('Hours');
subplot(3,1,[2])
plot(tHour, result_fixed(:,5), '.b:');
grid on; hold on;
ylim([0 max(abs(result_fixed(:,5)))]);
xlim([min(tHour) max(tHour)]);
ylabel('3D Error (meters)')
xlabel('Hours');
subplot(3,1,[3])
plot(tHour, result_fixed(:,8), '.b:');
grid on; hold on;
% ylim([0 0.05]);
xlim([min(tHour) max(tHour)]); 
ylabel('Number of Satellites')



