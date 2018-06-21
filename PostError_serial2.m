
dNE = rms(dNEV(:,4));
d3 = rms(dNEV(:,5));
tHour = [1:1:length(dNEV(:,1))]';

figure();
subplot(3,4,[1,2,5,6])
plot(dNEV(:,2), dNEV(:,1),'o'); 

axis([-5 5 -5 5]);

axis square
grid on; 
xlabel({'\Delta E (meters)',...
    ['dNE = ', num2str(dNE)],...
    [' 3D = ', num2str(d3)]}); ; ylabel('\Delta N (meters)')
%% 그래프 우측
subplot(3,4,[3,4])
plot(tHour(:,1), dNEV(:,1), '.r:', tHour(:,1), dNEV(:,2), '.b:'); xlim([1 length(dNEV)]); grid on;
legend('\Delta N', '\Delta E')
ylabel('\Delta N,E (meters)');
subplot(3,4,[7,8])
plot(tHour(:,1), dNEV(:,3), '.:'); xlim([1 length(dNEV)]); grid on;
ylabel('\Delta U (meters)')
