function [] = PlotVelDop(estm)


tth = gs2h24(estm(:,1));
tS = tth(1);
tE = tth(end);

vx = estm(:, 6)*3.6;
vy = estm(:, 7)*3.6;
vz = estm(:, 8)*3.6;
v3 = sqrt(vx.^2 + vy.^2 + vz.^2);

subplot(3,2,2); % Vx
plot(tth, vx, '.:')
xlim([tS tE]); ylabel('V_x');

subplot(3,2,4); % Vy
plot(tth, vy, '.:')
xlim([tS tE]); ylabel('V_y');

subplot(3,2,6); % Vz
plot(tth, vz, '.:')
xlim([tS tE]); ylabel('V_z'); xlabel('Hours')

subplot(3,2,3); % V
plot(tth, v3, '.:')
xlim([tS tE]); ylabel('V (km/h)');
xlabel('Hours');