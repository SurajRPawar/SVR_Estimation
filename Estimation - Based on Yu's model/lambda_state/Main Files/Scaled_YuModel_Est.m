% Simulates Yu's model, scaled model without estimation and a scaled model with estimation of the model.
% Script is written so that different versions of the estimation algorithm
% can be tried. For now, run version 5. Past tests also reveal promise with
% version 2
close all; clear; clc;

%% Parameters and Initial Conditions for Yu's simplified model
global Rc Lc Cao Rsvr Cr Ra Rm amp Rv;
vad_on = 0; % 1 = 5 L/min, 0 = 0 L/min
health = 1; % 1 = healthy, 0 = failure, Currently only tested and confirmed everything only for healthy case

parameters; % All parameters loaded from this script file.

% Select the version to simulate
% v1 : Estimating Rsvr, Cao, Cr using clubbed states
% v2 : Estimating Rsvr, Cao using clubbed states
% v3 : Estimating Rsvr, Cao, Cr, inverse states
% v4 : Estimating only Rsvr, Cr (under construction)
% v5 : Estimating only Rsvr 
% v6 : Estimating only Rsvr, but through inverse state (under construction)
version = 5;

%% Noise, Covariance parameters for estimation
Qk = 0.4; % Process noise, Higher Q -> Higher transients, less confidence on model
Pk = 0.2; % Initial error covariance. Higher P -> Lesser confidence in IC, slower transient

Q = diag(Qk.*[1,0.8,0.9,1.5]);
P = diag(Pk.*[1,0.001,1,0.001]);
R = (0.5)^2;


fprintf('Measurement uncertainty standard deviation +/- %.2f \n',sqrt(R));
fprintf('Measurement noise simulated is %.2f%% (?)\n',100*sqrt(R)/55);

% Q = diag(Qk.*[1,0.8,0.9,1.5]);


%% Scaling
% if vad_on == 0
%     Sx = diag([300, 50, 200, 2]);   % Scaling for states
% else
    Sx = diag([300, 50, 150, 1.4]);   % Scaling for states
% end

Sv0 = [5];                      % Scaling for initial volume V0
Sqvad = [5*1000/60];            % Scaling for QVAD measurements (treated as disturbance)

% Scale initial conditions
x0s = inv(Sx)*x0;
V0s = 1;
% P = inv(Sx)*P  % Scale P to match scaled state system
% Q = inv(Sx)*Q  % Scale Q to match scaled state system

fprintf('qao deviation %.2f%% \n',sqrt(P(1,1))*100/x0s(1));
fprintf('qr0 deviation %.2f%% \n',sqrt(P(2,2))*100/x0s(2));
fprintf('Vlv0 deviation %.2f%% \n',sqrt(P(3,3))*100/x0s(3));
fprintf('lambda0 deviation %.2f%% \n',sqrt(P(4,4))*100/x0s(4));
fprintf('Rsvr0 deviation %.2f%% \n \n',sqrt(Pk)*100/Rsvr);
%% Simulation parameters
HR   = 80;                    % [bpm]
Qvad_act = vad_on*5*1000/60;       % [mL/sec]
t0   = 0; tf = 60*10/HR;
dt   = 0.0002;                % For running model
tsamp = 1;              % For measurements
t    = [t0:dt:tf];
% tmeas = 
timesteps = length(t);

filt_cut_off = 5;   % [Hz]
filt_cut_off = filt_cut_off*2*pi;  %[rad/s]

filt_coeff = exp(-1*filt_cut_off*dt);

fprintf('Algorithm speed is %.3f kHz and measurements taken at %.3f kHz \n \n',1/(1000*dt),1/(1000*dt*tsamp));
%% Create Empty Vectors
% Unscaled system
x       = zeros(length(x0),timesteps); x(:,1) = x0;
zdp     = zeros(1,timesteps);
zplv    = zdp;
zpsa    = zdp;
e_t_lv  = zdp;

% Scaled system
xs = zeros(length(x0s),timesteps); xs(:,1) = x0s;
zdps    = zdp;
zplvs   = zdp;
zpsas   = zdp;
xis     = xs;   % Take scaled states, and get original. For comparison at the end of simulation.

% Estimation system
state0 = est_state(x0s,version, 1);
num_states = length(state0);
xhat    = zeros(num_states, timesteps); 

esthat = est_state(xhat(:,1),version,2);
num_ests = length(esthat);

% Append to noise matrices according to selected version
Q = blkdiag(Q,zeros(num_ests,num_ests));
P = blkdiag(P,Pk.*eye(num_ests));
P_all = zeros(length(P),length(P),timesteps);
P_all(:,:,1) = P;

% Randomize initialization of the initial states
rng('shuffle');   % Resets the random number generator so that results are different everytime
xhat(:,1) = state0 + chol(P)*randn(num_states,1);
zdphat  = zdp; zplvhat = zdp; zpsahat = zdp;
esthat = est_state(xhat(:,1),version,2);
num_ests = length(esthat);

Vtotalhat = xhat(1,1) + xhat(2,1) + xhat(3,1);

param_est = zeros(num_ests,timesteps); param_est(:,1) = esthat;

% Common vectors
DA_all = zeros(1,timesteps);    % Status of the mitral valve, 1 = Open, 0 = Closed
DM_all = zeros(1,timesteps);    % Status of the aortic valve, 1 = Open, 0 = Closed

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
Z_Qvads     = Z_Qvad./Sqvad;    % Scaled disturbance
fprintf('Interpolation of measurements took %.3f sec\n',toc)

% Add noise to measurements
dneu = chol(R)*randn(1,timesteps);
Z_deltaP_noisy = Z_deltaP + dneu;
Z_Qvads_noisy = Z_Qvads + dneu./Sqvad;

%% Parameters for Plv generation, taken from Jeff's paper
Emax = 3.25;                     % From Jeff's code
RR = 60/HR;                      % RR wave time [sec]
TLV = (1/1000)*(550 - 1.75*HR);  % Time for ventricular ejection [sec]

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

%% load ELV curve (from Ethan)
ELVn = load('From_CVS_Normalized_E.mat');
% ELVn  = load('ELV_normalized.mat');
time  = ELVn.time;
En    = ELVn.Enorm;
tmax0 = 0.2957;

Emax0 = 3.25; % Jeff's code
Emax = Emax0;
Ed0 = 0.006655; % Guess from Jeff's paper
Ed = Ed0; 
tc = 60/HR;
RR = tc;
ts = 0.07+0.3*tc; % systolic tirne interval % Ethan had this
% ts = 0.16+0.3*tc;
% ts = (1/1000)*(550 - 1.75*HR)  % Changing to this, based on Jeff's paper
TLV = ts;
% cardiovascular rmdel paranæters 

iter = 3; n = length(t);

tnvec = zeros(size(t));
tic
for k = 1:n
    t1 = rem(t(k),tc);
%     t1 = t(k);
    if t1 <= ts, tn = tmax0*(t1/ts); else; tn = tmax0 + (1-tmax0)*(t1-ts)/(tc-ts); end;
    tnvec(k) = tn;
%     Ev(k) = (Emax - Ed)*interp1(time,En,tn) + Ed;
end
Ev = (Emax - Ed)*interp1(time,En,tnvec) + Ed;
toc

%% Simulate Yu's model using Euler Integration
% Approach : Amplifying magnitude of valve resistances to simulate valve
% closure
% Jeff's code for Plv generation (active component only)
tic
for i = 1:timesteps   
    
    % Get elastance
    % E = Emax*e_t_lv(i);
    E = Ev(i);
    % Get current QVAD measurement
        Qvad = Z_Qvad(i);
    
    % Get xdot and outputs from function
    
        [fx, zs] = Yu(x(:,i),E,Qvad,V0);
    
    % Euler integration (expect at last step)
    if i ~= timesteps
        x(:,i+1) = x(:,i) + fx.*dt; 
    end
    
    % Collect outputs
        zdp(i)  = zs(1);
        zplv(i) = zs(2);
        zpsa(i) = zs(3); 
        DA_all(i) = zs(4);
        DM_all(i) = zs(5);
end
fprintf('Running unscaled model took %.3f sec \n',toc);

%% Simulate scaled Yu model
tic
for i = 1:timesteps
    
    % Get elastance
    % E = Emax*e_t_lv(i);
        E = Ev(i);
    % Get current QVAD measurement
        Qvad = Z_Qvads(i);
    
    % Get xdot and outputs from function
        [fx, zs] = Scaled_Yu(xs(:,i),E,Qvad,V0s);
    
    % Euler integration (expect at last step)
    if i ~= timesteps
        xs(:,i+1) = xs(:,i) + fx.*dt; 
    end
    
    % Collect outputs
        zdps(i)  = zs(1);
        zplvs(i) = zs(2);
        zpsas(i) = zs(3);
        DA_all(i) = zs(4);
        DM_all(i) = zs(5);
end
fprintf('Running scaled model took %.3f sec\n',toc);
xis = Sx*xs;
% makeplots_no_est;

%% Estimation with the scaled model
% Append based on number of estimates, which changes with the version used.
append = zeros(1,num_ests);
Rsvrhat_filt = zeros(1,timesteps);
Rsvrhat_filt(1) = param_est(1,1);

% In which iterations will measurements be available ?
trems = mod([1:timesteps], tsamp);
if tsamp ~=1, idx = trems == 1;
else idx = trems == 0;
end

Z_sparse = zeros(1,sum(idx));
Z_sparse(1) = Z_deltaP_noisy(1);
t_sparse = t(idx);
sparse_idx = 1;

fname = ['Yu_est_v',num2str(version)];
tic
for i = 2:timesteps
    
    %-------- Measurements at previous iteration -------------------------%
        Qvad = Z_Qvads_noisy(i-1);

    %-------- Measurements at current iteration  -------------------------%
        deltaP   = Z_deltaP_noisy(i);     
    
    %-------- State and error covariance propagation  --------------------%
        % Elastance at previous iteration
        E = Ev(i-1);
        % E = Emax*e_t_lv(i-1);
        % Model and covariance
        [fx_est, Pdot, zs] = feval(fname,xhat(:,i-1),E,Qvad,V0s,P,Q);
        
        xhat(:,i) = xhat(:,i-1) + fx_est.*dt;
        P = P + Pdot.*dt;
        
        if idx(i) == 1
            % Estimate only when aortic valve is closed (1)
            H = [0,0,0,7*Ra/(5*Lc),append]; 
            Qvad = Z_Qvads_noisy(i);
            KG = P*H.'/(H*P*H.' + R);
            P = (eye(num_states) - KG*H)*P;
            xhat(:,i) = xhat(:,i) + KG*(deltaP - 7*Ra*xhat(4,i)/(5*Lc) + 250*Qvad*Ra/3);
            Z_sparse(sparse_idx) = deltaP;
            sparse_idx = sparse_idx + 1;
        end
        
        % Adjust for total volume
        Vsum = xhat(1,i) + xhat(2,i) + xhat(3,i);
        %xhat(2,i) = xhat(2,i) + (Vtotalhat - Vsum);
        
        %----------- Gather all estimates --------------------------------%
        param_est(:,i) = est_state(xhat(:,i),version,2);
        
        %----------- Gather all outputs ----------------------------------%
        E = Ev(i);
        %E = Emax*e_t_lv(i);
        zdphat(i)  = 7*Ra*xhat(4,i)/(5*Lc) - 250*Z_Qvads_noisy(i)*Ra/3;
        zplvhat(i) = 150*E*xhat(3,i) - 5*E*V0s;
        zpsahat(i) = zplvhat(i) - zdphat(i);
        P_all(:,:,i) = P;
        
        %----------- Filtered Rsvr output --------------------------------%
        
        Rsvrhat_filt(i) = (1-filt_coeff)*param_est(1,i-1) + filt_coeff*Rsvrhat_filt(i-1);
end

fprintf('Running noisy model with estimation took %.3f sec\n',toc);

% Derive back original states from scaled system, to compare
xhat(1:4,:) = Sx*xhat(1:4,:);

%% Plots
makeplots;


%% Save data in file
% save('Scale_Yu.mat','t','x','zdp','zplv','zpsa','xs','zdps','zplvs','zpsas');
% fprintf('File saved successfully \n');