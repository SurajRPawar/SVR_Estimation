% Simulates Yu's model, scaled model without estimation and a scaled model with estimation.
% Script is written so that different versions of the estimation algorithm
% can be tried. Currently works with version 5.
close all; clear; clc;

%% Parameters and Initial Conditions for Yu's simplified model
global Rc Lc Cao Rsvr Cr Ra Rm amp Rv;
vad_on = 1; % 1 = 5 L/min, 0 = 0 L/min
health = 1; % 1 = healthy, 0 = failure

parameters; % Loads the electrical analog model's parameters, initial conditions (unscaled)

% Select the version to simulate
% v1 : Estimating Rsvr, Cao, Cr using clubbed states
% v2 : Estimating Rsvr, Cao using clubbed states
% v3 : Estimating Rsvr, Cao, Cr, inverse states
% v4 : Estimating only Rsvr, Cr (under construction)
% v5 : Estimating only Rsvr 
% v6 : Estimating only Rsvr, but through inverse state (under construction)
version = 5;

%% Noise, Covariance parameters for estimation
Qk = 40; % Process noise for scaled states, Higher Q -> Higher transients, less confidence on model
Pk = 1200; % Initial error covariance for scaled states. Higher P -> Lesser confidence in IC, slower transient

Q = diag(Qk.*[1,1,1,1]) % For the scaled states, non uniform distribution of the process noise to give more, or less weight to some states
P = diag(Pk.*[1.2,0.001,1,0.1]);    % This is for the scaled states
R = (5)^2;  % This is the uncertainty we want to place on measurements. Since I am using scaled states, tune this parameter heuristically.

%% Simulation parameters
HR   = 80;                    % [bpm]
Qvad_act = vad_on*5*1000/60;       % [mL/sec]
t0   = 0; tf = 60*200/HR;
dt   = 0.0002;                % For running model
tsamp = 10;              % For measurements
t    = [t0:dt:tf];
% tmeas = 
timesteps = length(t);

filt_cut_off = 1;   % [Hz]
filt_cut_off = filt_cut_off*2*pi;  %[rad/s]

filt_coeff = exp(-1*filt_cut_off*dt);

fprintf('Algorithm speed is %.3f kHz and measurements taken at %.3f kHz \n',1/(1000*dt),1/(1000*dt*tsamp));

%% Parameters for Plv generation, taken from Jeff's paper
Emax = 3.25;                     % From Jeff's code
RR = 60/HR;                      % RR wave time [sec]
% TLV = (1/1000)*(550 - 1.75*HR);  % Time for ventricular ejection [sec]

%% Elastance based on Jeff code
% Generate elastance curve
% for i = 1:timesteps
%            % Elastance
%        tn = mod(t(i),RR); 
%        if tn <= TLV*shift
%             e_t_lv(i) = 1/2 * (1 - cos(pi*tn/(TLV*shift)));
%         elseif tn > TLV*shift && tn <= TLV
%             e_t_lv(i) = 1/2 * (1 + cos(pi*(tn-TLV * shift)/(TLV*(1-shift))));
%         else
%             e_t_lv(i) = 0;
%        end
% end

%% load ELV curve (following Yu and Avonzollini's approach)
ELVn = load('From_CVS_Normalized_E.mat');
time  = ELVn.time;
En    = ELVn.Enorm;
tmax0 = 0.2957; % Time where Emax occurs. See on ELV vs t graph of one cycle with HR = 60bpm. 

Emax0 = 3.25; % Jeff's code
Emax = Emax0;
Ed0 = 0.006655; % Guess from Jeff's paper
Ed = Ed0; 
tc = 60/HR;
RR = tc;
ts = 0.07+0.3*tc; % systolic tirne interval, tuned in LabVIEW using Analysis_ELV_From_CVS.vi
% ts = 0.16+0.3*tc; % Avonzollini's equation mentioned in Yu's dissertation
% ts = (1/1000)*(550 - 1.75*HR)  % From Jeff's code
TLV = ts;

iter = 3; n = length(t);

tnvec = zeros(size(t));
tic
for k = 1:n
    t1 = rem(t(k),tc);
    if t1 <= ts, tn = tmax0*(t1/ts); else; tn = tmax0 + (1-tmax0)*(t1-ts)/(tc-ts); end;
    tnvec(k) = tn;
end
Ev = (Emax - Ed)*interp1(time,En,tnvec) + Ed;
toc

%% Load Measurement Data
savefile = ['CVS_12state_Zdata_',num2str(vad_on*5),'Qvad.mat'];
Zmat     = load(savefile);
Z_deltaP_all = Zmat.Z_deltaP;
Z_Qvad_all   = Zmat.Z_Qvad;
Z_PLV_all    = Zmat.Z_PLV;
Z_PSA_all    = Zmat.Z_PSA;
Z_time_all   = Zmat.teu;
Z_VLV_all = Zmat.Z_VLV;

% Interpolate measurements
% Create Measurement Vector
tic
Z_deltaP    = interp1(Z_time_all,Z_deltaP_all,t);
Z_Qvad      = interp1(Z_time_all,Z_Qvad_all,t);
Z_PLV       = interp1(Z_time_all,Z_PLV_all,t);
Z_VLV       = interp1(Z_time_all,Z_VLV_all,t);
Z_PSA       = interp1(Z_time_all,Z_PSA_all,t);
fprintf('Interpolation of measurements took %.3f sec\n',toc)

% Add noise to measurements
dneu = chol(R)*randn(1,timesteps);
% dneu = chol(R)*randn(1,timesteps);
Z_deltaP_noisy = Z_deltaP + dneu;
Z_Qvad_noisy = Z_Qvad + dneu;

%% Create Empty Vectors
% Unscaled system
x       = zeros(length(x0),timesteps); x(:,1) = x0;
zplv    = zeros(1,timesteps); zplv(1) = Ev(1)*(x0(3,1) - V0);
zpsa    = zplv; zpsa(1) = zplv(1) -Ra*(x0(4,1) - Z_Qvad(1));
zdp     = zplv; zdp(1) = zplv(1) - zpsa(1);

% Estimation system
state0 = est_state(x0,version, 1); % Populates the initial state for estimator (including estimates) without uncertainty
num_states = length(state0);
xhat    = zeros(num_states, timesteps); 

esthat = est_state(xhat(:,1),version,2);    % Extracts the first estimate from the initial estimate state
num_ests = length(esthat);

% Append to noise matrices according to selected version
Q = blkdiag(Q,zeros(num_ests,num_ests));
P = blkdiag(P,(Pk/10500).*eye(num_ests));
P_all = zeros(length(P),length(P),timesteps);
P_all(:,:,1) = P;

% Randomize initialization of the initial states
rng('shuffle');   % Resets the random number generator so that results are different everytime
xhat(:,1) = state0 + chol(P)*randn(num_states,1);   % Corrupt the initial estimate state by some uncertainty
zdphat  = zdp; zplvhat = zdp; zpsahat = zdp;
zplvhat(1) = Ev(1)*(xhat(3,1) - V0);
zpsahat(1) = zplvhat(1) - Ra*(xhat(4,1) - Z_Qvad_noisy(1));
zdphat(1) = zplvhat(1) - zpsahat(1);

esthat = est_state(xhat(:,1),version,2);    % Return back the uncertain initial parameter estimate

Vtotalhat = xhat(1,1) + xhat(2,1) + xhat(3,1);

param_est = zeros(num_ests,timesteps); param_est(:,1) = esthat;

% Common vectors
DA_all = zeros(1,timesteps);    % Status of the mitral valve, 1 = Open, 0 = Closed
DM_all = zeros(1,timesteps);    % Status of the aortic valve, 1 = Open, 0 = Closed

fprintf('qao std deviation %.2f%% \n',sqrt(P(1,1))*100/x0(1));
fprintf('qr0 std deviation %.2f%% \n',sqrt(P(2,2))*100/x0(2));
fprintf('Vlv0 std deviation %.2f%% \n',sqrt(P(3,3))*100/x0(3));
fprintf('fA0 std deviation %.2f%% \n',sqrt(P(4,4))*100/x0(4));
fprintf('Rsvr0 std deviation %.2f%% \n \n',sqrt(P(5,5))*100/Rsvr);

%% Simulate Yu's model using Euler Integration (clean measurements)
% Approach : Amplifying magnitude of valve resistances to simulate valve
% closure
tic
for i = 1:timesteps   
    
    % Get elastance
    % E = Emax*e_t_lv(i);  
    E = Ev(i);
    % Get current QVAD measurement
        Qvad = Z_Qvad(i);
    
    % Get xdot and outputs from function
        [fx, zs] = Yu_flowstate(x(:,i),E,Qvad,V0);
    
    % Euler integration (expect at last step)
    if i ~= timesteps
        x(:,i+1) = x(:,i) + fx.*dt; 
    end
    
    % Collect outputs
        zdp(i)  = zs(1);    % Differential pressure [mmHg]
        zplv(i) = zs(2);    % Pressure in LV [mmHg]
        zpsa(i) = zs(3);    % Pressure in systemic aorta [mmHg]
        DA_all(i) = zs(4);  % State of aortic valve, 1 = ON, Open; 0 = OFF, Closed
        DM_all(i) = zs(5);  % State of mitral valve, 1 = ON, Open; 0 = OFF, Closed
end

fprintf('Running unscaled model took %.3f sec\n',toc);

% makeplots_no_est;   % Plot comparing 12 state model, Yu model and Yu scaled model.

%% Estimation
% Append based on number of estimates, which changes with the version used.
append = zeros(1,num_ests); % This will be appended to the linearized H matrix
Rsvrhat_filt = zeros(1,timesteps);  % Empty array to store filtered estimates
Rsvrhat_filt(1) = param_est(1,1);   % Initialize the first filtered estimate to be same as unfiltered initial estimate

% In which iterations will measurements be available ?
trems = mod([1:timesteps], tsamp);  % Array of remainder when dividing the current iteration number, with the iteration at which measurement must be made available
if tsamp ~=1, idx = trems == 1;     % As long as you dont want algorithm and measurements to run at same rate, find the indexes (iterations) where remainder was 1.
else idx = trems == 0;              % If algorithm rate = measurement rate, look for indexes where remainder is 1.
end

Z_sparse = zeros(1,sum(idx));       % Create array to store the sparse measurements
Z_sparse(1) = Z_deltaP_noisy(1);    % Initialize first sparse measurement as the first measurement
t_sparse = t(idx);                  % Create the time vector containing time stamps when measurements will be provided to the estimator
sparse_idx = 1;                     % A counter to keep track of how many times a measurement update has been done.
KG_all = zeros(length(P),sum(idx));

fname = ['Yu_est_v',num2str(version),'_flowstate'];  % Based on defined version, change the function being called.
tic
for i = 2:timesteps
    
    %-------- Measurements at previous iteration -------------------------%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NOTE : When you are providing sparse measurements, you are still
    % restricted to get the Qvad measurement every iteration. Especially
    % true if the Qvad is pulsatile, or has transients. This makes it
    % important to have a Qvad state that your model can propagate as well.
    % For continuous flow pumps, you can probably provide this measurement
    % also sparsely. Not dealing with all that here. Assuming that I know
    % the VAD flow at all iteration steps.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Qvad = Z_Qvad_noisy(i-1);  

    %-------- Measurements at current iteration  -------------------------%
        deltaP   = Z_deltaP_noisy(i);     
    
    %-------- State and error covariance propagation  --------------------%
        % Elastance at previous iteration
        E = Ev(i-1);
        % E = Emax*e_t_lv(i-1);
        % Model and covariance
        [fx_est, Pdot, zs] = feval(fname,xhat(:,i-1),E,Qvad,V0,P,Q,DA_all(i),DM_all(i));
        
        xhat(:,i) = xhat(:,i-1) + fx_est.*dt;
        P = P + Pdot.*dt;
        
    %-------- Measurement update, only done when measurements available --%
        if idx(i) == 1
            H = [0,0,0,Ra,append]; % Scaled H matrix
%             H = [0,0,0,Ra/(Lc),append];
            Qvad = Z_Qvad_noisy(i);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % NOTE : In future, study how to simplify the calculations for
            % KG calculation, and P propagation. The equations as they are
            % now are difficult to implement on a controller, since they
            % involve huge matrix multiplications. Some ideas present in
            % the paper describing digital implementation of EKF for
            % motors.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            KG = P*H.'/(H*P*H.' + R);       % R, P all based on scaled states
            P = (eye(num_states) - KG*H)*P;
            xhat(:,i) = xhat(:,i) + KG*(deltaP - Ra*xhat(4,i) + Qvad*Ra);
            Z_sparse(sparse_idx) = deltaP;
            KG_all(:,sparse_idx) = KG;
            sparse_idx = sparse_idx + 1;
        end
        
        % Adjust for total volume
        Vsum = xhat(1,i) + xhat(2,i) + xhat(3,i);
        %xhat(2,i) = xhat(2,i) + (Vtotalhat - Vsum); % This step adjusts
        %the qr based on deviation from total volume. I thought this should
        %not be done, to allow the estimator to settle quickly. Let
        %estimator figure out where qr should be.
        
        %----------- Gather all estimates --------------------------------%
        param_est(:,i) = est_state(xhat(:,i),version,2);
        
        %----------- Gather all outputs ----------------------------------%
        E = Ev(i);
        %E = Emax*e_t_lv(i);
        zdphat(i)  = Ra*xhat(4,i) - Z_Qvad_noisy(i)*Ra;
        zplvhat(i) = E*xhat(3,i) - E*V0;
        zpsahat(i) = zplvhat(i) - zdphat(i);
        P_all(:,:,i) = P;
        
        %----------- Filtered Rsvr output --------------------------------%
        
        Rsvrhat_filt(i) = (1-filt_coeff)*param_est(1,i-1) + filt_coeff*Rsvrhat_filt(i-1);
end

fprintf('Running noisy model with estimation took %.3f sec\n',toc);

%% Plots
makeplots_flowstate;


%% Save data in file
% save('Scale_Yu.mat','t','x','zdp','zplv','zpsa','xs','zdps','zplvs','zpsas');
% fprintf('File saved successfully \n');