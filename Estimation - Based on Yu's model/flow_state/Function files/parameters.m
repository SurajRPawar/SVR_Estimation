%% CVS Parameters

% These are Yu's parameters
% amp =5e3;     % How much to amplify valve resistance to simulate valve closure
% Cao  = 1.33;       % [mL/mmHg]
% Rv  = 0.005;   % [mmHg.s/mL]
% Ra = 0.001;%Rv;        % [mmHg.s/mL]
% Rm = Rv;        % [mmHg.s/mL]
% Rsvr = 1;     % [mmHg.s/mL]    
% V0   = 5;       % Volume in LV at end of systole [mL]
% Cr = 4.4;       % [mL/mmHg]
% Rc = 0.0398;      % [mmHg.s/mL]
% Lc = 5e-4;    % [mmHg.s^2/mL]
% shift = 2/3;    % Used in normalized elastance calculation

% These are parameters tuned to match Jeff's impedance
amp =2e3;     % How much to amplify valve resistance to simulate valve closure
Cao  = 2;       % [mL/mmHg]
Rv  = 0.0025;   % [mmHg.s/mL]
Ra = Rv;        % [mmHg.s/mL]
Rm = Rv;        % [mmHg.s/mL]
Rsvr = 0.9;     % [mmHg.s/mL]    
V0   = 5;       % Volume in LV at end of systole [mL]
Cr = 0.8;       % [mL/mmHg]
Rc = 0.15;      % [mmHg.s/mL]
Lc = 2.2e-3;    % [mmHg.s^2/mL]
shift = 2/3;    % Used in normalized elastance calculation

%% Initial Conditions
% Values taken after running 12-state Gohean model
if vad_on == 1; %  VAD flow on at 5 L/min
    if health == 1;  % healthy case
        filename = 'IC_Chambers_5VAD_Healthy.mat';
    else
        filename = 'IC_Chambers_5VAD_Failure.mat';
    end
else
    if health == 1;  % healthy case
        filename = 'IC_Chambers_0VAD_Healthy.mat';
    else
        filename = 'IC_Chambers_0VAD_Failure.mat';
    end
end

ICs = load(filename);
PLV_0 = eval(['ICs.PLV0_',num2str(health)]);    % Left ventricle pressure [mmHg]
PLA_0 = eval(['ICs.PLA0_',num2str(health)]);     % Left atrium pressure [mmHg]
PRV_0 = eval(['ICs.PRV0_',num2str(health)]);     % Right ventricle pressure [mmHg]
PRA_0 = eval(['ICs.PRA0_',num2str(health)]);     % Right atrium pressure [mmHg]
PST_0 = eval(['ICs.PST0_',num2str(health)]);
VLV_0 = eval(['ICs.VLV0_',num2str(health)]);
QSA_0 = eval(['ICs.QSA0_',num2str(health)]);

% qao0 = Cao*(PST_0);   % Pst*Cao
% qr0 = (PLA_0)*Cr;
% Vlv0 = VLV_0;
% lambda0 = QSA_0*Lc;   % Qsa*Lc

if vad_on == 0;
    if health == 1;
        x0 = [159.0498; 6.1101; 115.5010; -14];%-0.0303];
    else
        x0 = [225.1180; 0.1500; 175.2530; QSA_0];%-0.0787];
    end
else
    if health == 1;
        x0 = [197.3748; 4.8576; 86.1588; 62.8718];%0.1402];
    else
        x0 = [274.2016; 0.1826; 159.4106; QSA_0];%0.0788];
    end
end

qao0    = x0(1);
qr0     = x0(2);
Vlv0    = x0(3);
fA0     = x0(4);
% lambda0 = x0(4);
Vtotal = Vlv0 + qr0 + qao0; 