% Makeplots_flowstate
%----------------- Figure 1 : Comparison of states -----------------------%
figure;
suptitle('States');
subplot(2,2,1);
plot(t,x(1,:),t,xhat(1,:),'--r'); grid on;
legend('Yu Orig','Estimated');
xlabel('Time [sec]'); ylabel('q_{ao} [mL]');
axisA(1) = gca;

subplot(2,2,2);
plot(t,Z_VLV,'k',t,x(3,:),'--b',t,xhat(3,:),'--r'); grid on;
xlabel('Time [sec]'); ylabel('V_{lv} [mL]');
legend('Jeff model','Yu Orig','Estimate');
axisA(2) = gca;

subplot(2,2,3);
plot(t,x(2,:),t,xhat(2,:),'--r'); grid on;
xlabel('Time [sec]'); ylabel('q_{r} [mL]');
legend('Yu Orig','Estimate');
axisA(3) = gca;

subplot(2,2,4);
plot(t,x(4,:),t,xhat(4,:),'--r'); grid on;
xlabel('Time [sec]'); ylabel('f_A [mL/sec]');
legend('Yu Orig','Estimate');
axisA(4) = gca;
linkaxes(axisA,'x');

%-------------------- Figure 2 : Comparison of measurements --------------%
figure; 
suptitle('Outputs');
subplot(4,1,1)
plot(t,Z_PLV,t,zplv,'--b',t,zplvhat,'.','MarkerSize',6); grid on;
xlabel('Time [sec]'); ylabel('P_{lv} [mmHg]');
title('P_{lv} comparison');
legend('Jeff model','Yu Orig','Estimate');
axisB(1) = gca;

subplot(4,1,2)
plot(t,Z_PSA,t,zpsa,'--b',t,zpsahat,'.','MarkerSize',6); grid on;%,t,zpsas,'--r'
xlabel('Time [sec]'); ylabel('P_{sa} [mmHg]');
title('P_{sa} comparison');
legend('Jeff model','Yu Orig','Estimate');
axisB(2) = gca;

subplot(4,1,3)
plot(t,Z_deltaP,t,zdp,'--b',t,zdphat,'.','MarkerSize',6); grid on;
xlabel('Time [sec]'); ylabel('\Delta_P [mmHg]');
title('Delta Pressure measurement comparison');
legend('Jeff model','Yu Orig','Estimate');
axisB(3) = gca;

subplot(4,1,4)
plot(t_sparse,Z_sparse,'o'); grid on;
xlabel('Time [sec]'); ylabel('\Delta_P [mmHg]');
title('Sparse Delta Pressure measurement');
axisB(4) = gca;

linkaxes(axisB,'x');

if version == 1 || version == 3
    Rsvrhat = param_est(1,:);
    Caohat = param_est(2,:);
    Crhat = param_est(3,:);
    %-------------------- Figure 3 : EKF parameter estimates -------------%
    figure;
    suptitle('Parameter Estimates');
    
    subplot(3,1,1);
    plot(t,Rsvr*ones(1,timesteps),t,Rsvrhat,'--r'); grid on;
    xlabel('Time [sec]'); ylabel('R_{svr} [mmHg.s/mL]');
    legend('Actual','Estimate');
    
    subplot(3,1,2);
    plot(t,Cao*ones(1,timesteps),t,Caohat,'--r'); grid on;
    xlabel('Time [sec]'); ylabel('C_{ao} [mL/mmHg]');
    legend('Actual','Estimate');
    
    subplot(3,1,3);
    plot(t,Cr*ones(1,timesteps),t,Crhat,'--r'); grid on;
    xlabel('Time [sec]'); ylabel('C_{r} [mL/mmHg]');
    legend('Actual','Estimate');
    
elseif version == 2
    Rsvrhat = param_est(1,:);
    Caohat = param_est(2,:);
    
    figure;
    suptitle('Parameter Estimates');
    
    subplot(2,1,1);
    plot(t,Rsvr*ones(1,timesteps),t,Rsvrhat,'--r'); grid on;
    xlabel('Time [sec]'); ylabel('R_{svr} [mmHg.s/mL]');
    legend('Actual','Estimate');
    
    subplot(2,1,2);
    plot(t,Cao*ones(1,timesteps),t,Caohat,'--r'); grid on;
    xlabel('Time [sec]'); ylabel('C_{ao} [mL/mmHg]');
    legend('Actual','Estimate');
elseif version == 5
    Rsvrhat = param_est(1,:);
    Rsvrhat_u = 1.1*Rsvr;
    Rsvrhat_l = 0.9*Rsvr;
    figure;
    subplot(2,1,1);
    plot(t,Rsvr*ones(1,timesteps),t,Rsvrhat,'--r',t,Rsvrhat_u*ones(1,timesteps),'-b',t,Rsvrhat_l*ones(1,timesteps),'-b'); grid on;
    xlabel('Time [sec]'); ylabel('R_{svr} [mmHg.s/mL]');
    legend('Actual','Estimate');
    title('R_{svr} estimate');
    
    subplot(2,1,2);
    plot(t,Rsvrhat_filt); grid on;
    xlabel('Time [sec]'); ylabel('R_{svr} [mmHg.s/mL]');
    legend('Filtered');
    title('R_{svr} estimate');
    
    % ----------------------- Figure 4 : Error Covariance-----------------%
    P_qao = squeeze(P_all(1,1,:));
    P_qr = squeeze(P_all(2,2,:));
    P_Vlv = squeeze(P_all(3,3,:));
    P_fA = squeeze(P_all(4,4,:));
    P_Rsvr = squeeze(P_all(5,5,:));

    figure;
    suptitle(['Error Covariance at T_s =',num2str(dt*1000),'ms']);
    subplot(3,2,1);
    plot(t,P_qao); 
    xlabel('Time [sec]'); ylabel('Error Covariance in q_ao');

    subplot(3,2,2);
    plot(t,P_qr); 
    xlabel('Time [sec]'); ylabel('Error Covariance in q_r');

    subplot(3,2,3);
    plot(t,P_Vlv); 
    xlabel('Time [sec]'); ylabel('Error Covariance in V_{lv}');

    subplot(3,2,4);
    plot(t,P_fA); 
    xlabel('Time [sec]'); ylabel('Error Covariance in f_A');

    subplot(3,2,[5,6]);
    plot(t,P_Rsvr); 
    xlabel('Time [sec]'); ylabel('Error Covariance in R_{svr}');
    
    % ----------------------- Figure 5 : Kalman Gains    -----------------%
    KG_qao = KG_all(1,:);
    KG_qr =  KG_all(2,:);
    KG_Vlv = KG_all(3,:);
    KG_fA =  KG_all(4,:);
    KG_Rsvr = KG_all(5,:);

    figure;
    suptitle(['Kalman Gains at T_s =',num2str(dt*1000),'ms']);
    subplot(3,2,1);
    plot(t_sparse,KG_qao); 
    xlabel('Time [sec]'); ylabel('Kalman Gain q_ao');

    subplot(3,2,2);
    plot(t_sparse,KG_qr); 
    xlabel('Time [sec]'); ylabel('Kalman Gain q_r');

    subplot(3,2,3);
    plot(t_sparse,KG_Vlv); 
    xlabel('Time [sec]'); ylabel('Kalman Gain V_{lv}');

    subplot(3,2,4);
    plot(t_sparse,KG_fA); 
    xlabel('Time [sec]'); ylabel('Kalman Gain f_A');

    subplot(3,2,[5,6]);
    plot(t_sparse,KG_Rsvr); 
    xlabel('Time [sec]'); ylabel('Kalman Gain R_{svr}');
end