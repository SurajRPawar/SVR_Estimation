% Make plots just to compare Jeff's model, Yu model and unscaled model.
figure; 
subplot(4,1,1)
plot(t,Z_PLV,t,zplv,'--b',t,zplvs,'--r'); grid on;
xlabel('Time [sec]'); ylabel('P_{lv} [mmHg]');
title('P_{lv} comparison');
legend('Jeff model','Yu Orig','Yu Scaled');
axisB(1) = gca;

subplot(4,1,2)
plot(t,Z_PSA,t,zpsa,'--b',t,zpsas,'--r'); grid on;%,t,zpsas,'--r'
xlabel('Time [sec]'); ylabel('P_{sa} [mmHg]');
title('P_{sa} comparison');
legend('Jeff model','Yu Orig','Yu Scaled');%,'Yu Scaled'
axisB(2) = gca;

subplot(4,1,3)
plot(t,Z_deltaP,t,zdp,'--b',t,zdps,'--r'); grid on;
xlabel('Time [sec]'); ylabel('\Delta_P [mmHg]');
title('Delta Pressure measurement comparison');
legend('Jeff model','Yu Orig','Yu Scaled');
axisB(3) = gca;

subplot(4,1,4);
plot(t,Z_VLV,'k',t,x(3,:),'--b',t,xis(3,:),'--r'); grid on;
xlabel('Time [sec]'); ylabel('V_{lv} [mL]');
legend('Jeff model','Yu Orig','Yu Scaled');
axisB(4) = gca;

linkaxes(axisB,'x');