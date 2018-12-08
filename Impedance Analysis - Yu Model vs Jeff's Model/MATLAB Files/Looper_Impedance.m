close all; clear all; clc;

global Ca
Cao  = 1.33;       % [mL/mmHg]
Rv  = 0.0005;   % [mmHg.s/mL]
Ra = Rv;        % [mmHg.s/mL]
Rm = Rv;        % [mmHg.s/mL]
Rsvr = 0.9;       % [mmHg.s/mL]    
Vo   = 5;       % Volume in LV at end of systole [mL]
Cr = 4.4;         % [mL/mmHg]
Rc = 0.0398;      % [mmHg.s/mL]
Lc = 5e-4;%2.2e-3;    % [mmHg.s^2/mL]
rate_hz = [1:0.5:5];
rate = rate_hz.*60;
steps = length(rate);
for i = 1:steps
    MCL_CVS_Gohean_function_version(rate(i));
end