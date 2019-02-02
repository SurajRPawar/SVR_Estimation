%% Test MCL_CVS_Gohean_ODE.m
close all; clear all; clc

%% Variables
tend = 60/80*100;
dt   = 0.0001;
t    = [0:dt:tend];
rho  = 1.055*1000; %[kg/m^3]
mu   = 0.0027;
HR   = 80;
Qvad = 0*1000/60; % [L/min]*[1000mL/L]*[min/60s] = [mL/s]

% CVS states: x = [PSA;PST;PSV;PPA;PPT;PPV;QSA;QPA;VLA;VLV;VRA;VRV];
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
    5.4813]

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
   18.0235]

% x0cvs = [61.5086  118.7781   86.2430   38.0359   80.5853    6.3602   43.4843  110.0420 14.2956   64.2535    9.7965    6.7447];
% x0cvs = [x0cvs(3)
%     x0cvs(5)
%     x0cvs(6)
%     x0cvs(9)
%     x0cvs(11)
%     x0cvs(12)
%     x0cvs(4)
%     x0cvs(10)
%     x0cvs(1)
%     x0cvs(2)
%     x0cvs(7)
%     x0cvs(8)]

x0 = [x0_c1, x0_c2];

%% Cases: 1 = Healthy; 2 = CHF
cases  = 2; % Total number of cases
Case   = [1 2]; % Chosen case #

%% Simulation
options = odeset('MaxStep',0.001);
TmL     = zeros(length(t),cases); XmL = zeros(length(TmL),length(x0),cases);
t_short = [0:dt:60/80*20]; % 20 cycles to get to steady state

tic
% Find Steady State CVS ICs
for i = 1:cases
    [T,X]   = ode23(@MCL_CVS_Gohean_ODE,t_short,x0(:,i),options,rho,mu,HR,Case(i),Qvad); 
    X_short = X;
    
    % Update ICs
    x0_new(:,i) = X_short(end,:)'
end

% Run Simulation with New ICs
for i = 1:cases
    [T,X]      = ode23(@MCL_CVS_Gohean_ODE,t,x0_new(:,i),options,rho,mu,HR,Case(i),Qvad);
    TmL(:,i)   = T; 
    XmL(:,:,i) = X;
end
toc

%% Outputs
VLVmL1 = XmL(:,2,1); VLVmL2 = XmL(:,2,2);

PSA = zeros(size(TmL));
PLV = PSA;
% VLAd = PSA; VLVd = PSA; VSAd = PSA; VATd = PSA;

tic
it0 = 1;
for j = 1:cases
    for i = it0:length(TmL(:,j))
        [xdot y]  = MCL_CVS_Gohean_ODE(TmL(i,j),XmL(i,:,j),rho,mu,HR,Case(j),Qvad);
        PSA(i,j)  = y(1);
        PLV(i,j)  = y(2);
    end
end
toc

PLVmL1 = PLV(:,1); PLVmL2 = PLV(:,2);
PSAmL1 = PSA(:,1); PSAmL2 = PSA(:,2);

%% Euler
% dt    = 0.001;
teu   = [0:dt:tend];
Xeu   = zeros(length(teu),12,cases); Xeu(1,:,1) = x0_new(:,1)'; Xeu(1,:,2) = x0_new(:,2)';
PLVeu = zeros(length(teu),cases); PSAeu = PLVeu; ELV = PSAeu; Qsv = ELV;
tic

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
    end
end
toc

PSAeu1 = PSAeu(:,1); PSAeu2 = PSAeu(:,2);
PLVeu1 = PLVeu(:,1); PLVeu2 = PLVeu(:,2);
VLVeu1 = Xeu(:,2,1); VLVeu2 = Xeu(:,2,2);

%% Figures
figure;
plot(t,PLVmL1,t,PSAmL1,teu,PLVeu1,'--',teu,PSAeu1,'--'); grid on;
xlabel('Time [s]'); ylabel('Pressure [mmHg]');
legend('PLV','PSA','PLVeu','PSAeu');
title(['Pressure generation with ode23 and euler; qVAD = ',num2str(Qvad*60/1000),' [LPM]']);

figure;
i1 = 10001; i2 = 15001; 
plot(VLVmL1,PLVmL1,VLVeu1,PLVeu1,'--',VLVmL2,PLVmL2,'k',VLVeu2,PLVeu2,'--'); grid on
% plot(VLVmL1(i1:i2),PLVmL1(i1:i2),VLVeu1(i1:i2),PLVeu1(i1:i2),'--',VLVmL2(i1:i2),PLVmL2(i1:i2),'k',VLVeu2(i1:i2),PLVeu2(i1:i2),'--'); grid on
axis([0 300 0 150])
xlabel('LV Volume [mL]'); ylabel('LV Pressure [mmHg]');
legend('healthy ode23','healthy euler','failure ode23','failure euler')
title(['P-V Loops with ode23 and euler; qVAD = ',num2str(Qvad*60/1000),' [LPM]']);

%% DeltaP for Measurement Data (taken after 10 cycles)
Tend = 80*60*10;
Z_deltaP = (PLVeu1 - PSAeu1)';
Z_Qvad   = Qvad*ones(size(teu));
Z_PLV    = PLVeu1';
Z_PSA    = PSAeu1';
figure;
subplot(2,1,1)
plot(teu,Z_deltaP); grid on;
ylabel('\Delta_P [mmHg]')

subplot(2,1,2)
plot(teu,Z_Qvad); grid on;
ylabel('Q_{VAD} [L/s]'); xlabel('Time [s]');

%% Save Measurement Data
% SaveData = 'No';
SaveData = 'Yes';

TF = strcmp(SaveData,'Yes');
if TF == 1;
    savefile = ['CVS_12state_Zdata_',num2str(Qvad*60/1000),'Qvad.mat'];
    save(savefile,'teu','Z_deltaP','Z_Qvad','Z_PLV','Z_PSA');
    fprintf('\n Measurement Data Saved\n');
    % Check Saved File Data
    Zmat = load(savefile);
    Z_deltaP = Zmat.Z_deltaP;
    Z_Qvad   = Zmat.Z_Qvad;
    Z_PLV    = Zmat.Z_PLV;
    Z_PSA    = Zmat.Z_PSA;
else
    fprintf('\n Measurement Data Not Saved\n');
end






















