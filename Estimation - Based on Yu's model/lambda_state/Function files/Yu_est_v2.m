% Version 2 - Estimating Rsvr, Cao. Cr assumed known
function [fx_est, Pdot, y] = Yu_est_v1(xhat,E,Qvad,V0s,P,Q);
global Rc Lc Cao Rsvr Cr Ra Rm amp Rv;

% States and estimates
qaohat    = xhat(1);
qrhat     = xhat(2);
Vlvhat    = xhat(3);
lambdahat = xhat(4);
x5hat     = xhat(5);
x6hat     = xhat(6);

Rsvrhat = 1/(Cr*x5hat);
Caohat = 1/x6hat;

% Variables to calculate
fA = 2*lambdahat/Lc;
Plv = E*(200*Vlvhat - 5*V0s);
Psa = 250*Qvad*Ra/3 + Plv - Ra*fA;
Pla = 50*qrhat/Cr; 

% Toggle valve status
if Plv > Psa
    % Ventricle Ejection stage, Aortic valve open (1)
    DA = 1;
    Ra = Rv;
else
    % Ventricle filling stage, Aortic valve closed (0)
    DA = 0;
    Ra = amp*Rv;
end

if Pla > Plv
    % Atria ejection stage, Mitral valve open (1)
    DM = 1;
    Rm = Rv;
else
    DM = 0;
    Rm = amp*Rv;
end

fx_est = [lambdahat/(150*Lc) - qaohat*x5hat*x6hat*Cr + qrhat*x5hat/6; 
          6*qaohat*x5hat*x6hat*Cr - qrhat/(Cr*Rm) - qrhat*x5hat + 4*E*Vlvhat/Rm - E*V0s/(10*Rm);
          qrhat/(4*Cr*Rm) - lambdahat/(100*Lc) - E*Vlvhat/Rm + E*V0s/(40*Rm);
          125*Qvad*Ra/3 - 5*E*V0s/2 - 150*qaohat*x6hat + 100*E*Vlvhat - lambdahat*(Ra+Rc)/Lc;
          0;
          0];
      
F = ...
[  -Cr*x5hat*x6hat,             x5hat/6,        0,    1/(150*Lc), qrhat/6 - Cr*qaohat*x6hat,  -Cr*qaohat*x5hat;
  6*Cr*x5hat*x6hat, - x5hat - 1/(Cr*Rm), (4*E)/Rm,             0, 6*Cr*qaohat*x6hat - qrhat, 6*Cr*qaohat*x5hat;
                 0,         1/(4*Cr*Rm),    -E/Rm,   -1/(100*Lc),                         0,                 0;
        -150*x6hat,                   0,    100*E, -(Ra + Rc)/Lc,                         0,       -150*qaohat;
                 0,                   0,        0,             0,                         0,                 0;
                 0,                   0,        0,             0,                         0,                 0];
      
Pdot = F*P + P*F.' + Q;

% Outputs
y(1) = Plv - Psa; % Delta P [mmHg]
y(2) = Plv; % Left Ventricle pressure [mmHg]
y(3) = Psa; % Systemic Aorta pressure [mmHg]
y(4) = DA;  % Aortial valve status
y(5) = DM;  % Mitral valve status
end