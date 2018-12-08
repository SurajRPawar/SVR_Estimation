function MCL_CVS_Gohean_function_version(Rate)
% CVS states: x = VLA;VLV;P_sa;Q_sa;P_st;P_sv;VRA;VRV;P_pa;Q_pa;P_pt;P_pv;
% x0_new =
% 
%    29.7433   73.8699
%    88.4510  264.9593
%   113.9210   95.0877
%    96.0422   83.2879
%    99.5350   82.5972
%     8.4634    7.8827
%    52.9784   54.6765
%   115.1204  116.6454
%    11.1595   18.0337
%    81.9227   60.8742
%     5.5552   13.7383
%     1.9785   10.9797

%% Test MCL_CVS_Gohean_ODE.m
% close all; clear all; clc

%% Variables
HR   = Rate;
tend = 60/HR*50;
dt   = 0.0001;
t    = [0:dt:tend];
rho  = 1.055*1000; %[kg/m^3]
mu   = 0.0027;
Qvad = 0*5*1000/60; % [L/min]*[1000mL/L]*[min/60s] = [mL/s]

% CVS states: x = VLA;VLV;P_sa;Q_sa;P_st;P_sv;VRA;VRV;P_pa;Q_pa;P_pt;P_pv;
x0_c1 = [55.5739
  115.3849
   87.9005
   39.3922
   82.0488
    5.8710
   38.9063
  106.1946
   13.5506
   69.1017
    8.7373
    5.4813];

x0_c2 = [84.2733
  293.8926
   48.4854
   14.9607
   46.2601
    3.3689
   30.5154
  100.4741
   22.3319
   39.6700
   19.6512
   18.0235];

x0 = [x0_c1, x0_c2];

%% Cases: 1 = Healthy; 2 = CHF
cases  = 2; % Total number of cases
Case   = [1 2]; % Chosen case #

%% Simulation
options = odeset('MaxStep',0.001);
TmL     = zeros(length(t),cases); XmL = zeros(length(TmL),length(x0),cases);
t_short = [0:dt:60/80*20]; % 20 cycles to get to steady state

% tic
% Find Steady State CVS ICs
for i = 1:cases
    [T,X]   = ode23(@MCL_CVS_Gohean_ODE,t_short,x0(:,i),options,rho,mu,HR,Case(i),Qvad); 
    X_short = X;
    
    % Update ICs
    x0_new(:,i) = X_short(end,:)';
end

% Run Simulation with New ICs
% for i = 1:cases
%     [T,X]      = ode23(@MCL_CVS_Gohean_ODE,t,x0_new(:,i),options,rho,mu,HR,Case(i),Qvad);
%     TmL(:,i)   = T; 
%     XmL(:,:,i) = X;
% end
% toc

%% Outputs
VLVmL1 = XmL(:,2,1); VLVmL2 = XmL(:,2,2);

PSA = zeros(size(TmL));
PLV = PSA;
PLA = PSA;
PRV = PSA;
PRA = PSA;

% VLAd = PSA; VLVd = PSA; VSAd = PSA; VATd = PSA;

% tic
it0 = 1;
for j = 1:cases
    for i = it0:length(TmL(:,j))
        [xdot y]    = MCL_CVS_Gohean_ODE(TmL(i,j),XmL(i,:,j),rho,mu,HR,Case(j),Qvad);
        PSA(i,j)    = y(1);
        PLV(i,j)    = y(2);
        e_t_lv(i,j) = y(3);
        PLA(i,j) = y(4);
        PRV(i,j) = y(8);
        PRA(i,j) = y(9);
    end
end
% toc

PLVmL1 = PLV(:,1); PLVmL2 = PLV(:,2);
PSAmL1 = PSA(:,1); PSAmL2 = PSA(:,2);
PLAmL1 = PLA(:,1); PLAmL2 = PLA(:,2);
PRAmL1 = PRA(:,1); PRAmL2 = PRA(:,2);
PRVmL1 = PRV(:,1); PRVmL2 = PRV(:,2);

%% Euler
% dt    = 0.001;
teu   = [0:dt:tend];
Xeu   = zeros(length(teu),12,cases); Xeu(1,:,1) = x0_new(:,1)'; Xeu(1,:,2) = x0_new(:,2)';
PLVeu = zeros(length(teu),cases); PSAeu = PLVeu; ELV = PSAeu; Qsv = ELV; PLAeu = PSAeu;
PRVaeu = PLVeu;
PRAaeu = PLVeu;
PLVa = PLVeu;

% tic

for j = 1:cases
    for i = 1:length(teu)
        if i == 1
            [xdot,y] = MCL_CVS_Gohean_ODE(teu(i),Xeu(i,:,j)',rho,mu,HR,Case(j),Qvad);
        else
            [xdot,y] = MCL_CVS_Gohean_ODE(teu(i-1),Xeu(i-1,:,j)',rho,mu,HR,Case(j),Qvad);
            Xeu(i,:,j) = Xeu(i-1,:,j) + xdot'.*dt;
        end
        PSAeu(i,j) = y(1);
        PLVeu(i,j) = y(2);
        PLAeu(i,j) = y(4);
        PLAa(i,j) = y(5);
        PRVa(i,j) = y(6);
        PRVeu(i,j) = y(8);
        PRAa(i,j) = y(7);
        PRAeu(i,j) = y(9);
        PLVa(i,j) = y(10);
    end
end
% toc

PSAeu1 = PSAeu(:,1); PSAeu2 = PSAeu(:,2);
PLVeu1 = PLVeu(:,1); PLVeu2 = PLVeu(:,2);
VLVeu1 = Xeu(:,2,1); VLVeu2 = Xeu(:,2,2);
PLAeu1 = PLAeu(:,1); PLAeu2 = PLAeu(:,2);
PRVeu1 = PRVeu(:,1); PRVeu2 = PRVeu(:,2);
PRAeu1 = PRAeu(:,1); PRAeu2 = PRAeu(:,2);

%% Figures
% figure;
% plot(t,PLVmL1,t,PSAmL1,teu,PLVeu1,'--',teu,PSAeu1,'--'); grid on;
% xlabel('Time [s]'); ylabel('Pressure [mmHg]');
% legend('PLV','PSA','PLVeu','PSAeu');
% title(['Pressure generation with ode23 and euler; qVAD = ',num2str(Qvad*60/1000),' [LPM]']);
% 
% figure;
% i1 = 10001; i2 = 15001; 
% plot(VLVmL1,PLVmL1,VLVeu1,PLVeu1,'--',VLVmL2,PLVmL2,'k',VLVeu2,PLVeu2,'--'); grid on
% % plot(VLVmL1(i1:i2),PLVmL1(i1:i2),VLVeu1(i1:i2),PLVeu1(i1:i2),'--',VLVmL2(i1:i2),PLVmL2(i1:i2),'k',VLVeu2(i1:i2),PLVeu2(i1:i2),'--'); grid on
% axis([0 300 0 150])
% xlabel('LV Volume [mL]'); ylabel('LV Pressure [mmHg]');
% legend('healthy ode23','healthy euler','failure ode23','failure euler')
% title(['P-V Loops with ode23 and euler; qVAD = ',num2str(Qvad*60/1000),' [LPM]']);
% 
% %% DeltaP for Measurement Data (taken after 10 cycles)
% Tend = 80*60*10;
% Z_deltaP = (PLVeu1 - PSAeu1)';
% Z_Qvad   = Qvad*ones(size(teu));
% Z_PLV    = PLVeu1';
% Z_PSA    = PSAeu1';
% Z_PLA = PLAeu1.';
% Z_PRV = PRVeu(:,1).';
% Z_PRA = PRAeu(:,1).';
% Z_VLV = Xeu(:,2,1).';
% 
% figure;
% subplot(2,1,1)
% plot(teu,Z_deltaP); grid on;
% ylabel('\Delta_P [mmHg]')
% 
% subplot(2,1,2)
% plot(teu,Z_Qvad); grid on;
% ylabel('Q_{VAD} [L/s]'); xlabel('Time [s]');
% 
% %% Save Measurement Data
% SaveData = 'Yes';
% % SaveData = 'Yes';
% 
% TF = strcmp(SaveData,'Yes');
% if TF == 1;
%     savefile = ['CVS_12state_Zdata_',num2str(Qvad*60/1000),'Qvad.mat'];
%     save(savefile,'teu','Z_deltaP','Z_Qvad','Z_PLV','Z_PSA','Z_PLA', 'Z_PRV','Z_PRA','Z_VLV');
%     fprintf('\n Measurement Data Saved\n');
%     % Check Saved File Data
%     Zmat = load(savefile);
%     Z_deltaP = Zmat.Z_deltaP;
%     Z_Qvad   = Zmat.Z_Qvad;
%     Z_PLV    = Zmat.Z_PLV;
%     Z_PSA    = Zmat.Z_PSA;
% else
%     fprintf('\n Measurement Data Not Saved\n');
% end
% 
% %% Elastance Figures
% ELVmax = 3.25; V0 = 5; 
% figure;
% % subplot(2,1,1)
% plot(t,e_t_lv.*ELVmax); grid on;
% 
% figure;
% % subplot(2,1,2)
% 
% plot(t,(e_t_lv.*ELVmax.*(VLVeu1 - V0)),t,PLVeu1,'--r'); grid on;
% legend('Passive PLV','Actual PLV');
% 
% %% Comparison of active and total pressures - RV, RA, LA
% figure;
% subplot(2,2,1);
% plot(t,PLVeu1,t,PLVa(:,1),'--r');
% xlabel('Time [sec]'); ylabel('Pressure [mmHg]');
% legend('Total','Active');
% title('Left Ventricle Pressure');
% 
% subplot(2,2,2);
% plot(t,PLAeu1,t,PLAa(:,1),'--r');
% xlabel('Time [sec]'); ylabel('Pressure [mmHg]');
% legend('Total','Active');
% title('Left Atrium Pressure');
% 
% subplot(2,2,3);
% plot(t,PRVeu(:,1),t,PRVa(:,1),'--r');
% xlabel('Time [sec]'); ylabel('Pressure [mmHg]');
% legend('Total','Active');
% title('Right Ventricle Pressure');
% 
% subplot(2,2,4);
% plot(t,PRAeu(:,1),t,PRAa(:,1),'--r');
% xlabel('Time [sec]'); ylabel('Pressure [mmHg]');
% legend('Total','Active');
% title('Right Atrium Pressure');

% figure;
% subplot(2,1,1);
% plot(teu,Xeu(:,4,1),'.');
% xlabel('Time [sec]'); ylabel('Q_{sa} [mL/sec]');
% title('Flow in artery');
% 
% subplot(2,1,2);
% plot(teu,Xeu(:,3,1),'.');
% xlabel('Time [sec]'); ylabel('P_{sa} [mmHg]');
% title('Pressure in artery');

filename = ['Impedance_Jeff_',num2str(HR),'bpm.txt'];
M = [teu', Xeu(:,3,1), Xeu(:,4,1)];
csvwrite(filename,M);
end















