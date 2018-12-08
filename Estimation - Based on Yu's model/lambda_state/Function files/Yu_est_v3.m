% Version 3 - Estimating inverses of parameters
function [fx_est, Pdot, y] = Yu_est_v1(xhat,E,Qvad,V0s,P,Q);
global Rc Lc Cao Rsvr Cr Ra Rm amp Rv;

% States and estimates
qaohat    = xhat(1);
qrhat     = xhat(2);
Vlvhat    = xhat(3);
lambdahat = xhat(4);
x5hat     = xhat(5);
x6hat     = xhat(6);
x7hat     = xhat(7);

Rsvrhat = 1/x5hat;
Caohat  = 1/x6hat;
Crhat   = 1/x7hat;

% Variables to calculate
fA = 2*lambdahat/Lc;
Plv = E*(200*Vlvhat - 5*V0s);
Psa = 250*Qvad*Ra/3 + Plv - Ra*fA;
Pla = 50*qrhat/Crhat; 

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

fx_est = [lambdahat/(150*Lc) - qaohat*x5hat*x6hat + qrhat*x5hat*x7hat/6; 
          6*qaohat*x5hat*x6hat - qrhat*(x7hat/Rm + x5hat*x7hat) + 4*E*Vlvhat/Rm - E*V0s/(10*Rm);
          qrhat*x7hat/(4*Rm) - lambdahat/(100*Lc) - E*Vlvhat/Rm + E*V0s/(40*Rm);
          125*Qvad*Ra/3 - 5*E*V0s/2 - 150*x6hat*qaohat + 100*E*Vlvhat - lambdahat*(Ra+Rc)/Lc;
          0;
          0;
          0];
      
F = ...
[  -x5hat*x6hat,          (x5hat*x7hat)/6,        0,    1/(150*Lc), (qrhat*x7hat)/6 - qaohat*x6hat,  -qaohat*x5hat,       (qrhat*x5hat)/6;
  6*x5hat*x6hat, - x5hat*x7hat - x7hat/Rm, (4*E)/Rm,             0,   6*qaohat*x6hat - qrhat*x7hat, 6*qaohat*x5hat, -qrhat*(x5hat + 1/Rm);
              0,             x7hat/(4*Rm),    -E/Rm,   -1/(100*Lc),                              0,              0,          qrhat/(4*Rm);
     -150*x6hat,                        0,    100*E, -(Ra + Rc)/Lc,                              0,    -150*qaohat,                     0;
              0,                        0,        0,             0,                              0,              0,                     0;
              0,                        0,        0,             0,                              0,              0,                     0;
              0,                        0,        0,             0,                              0,              0,                     0];
      
Pdot = F*P + P*F.' + Q;

% Outputs
y(1) = Plv - Psa; % Delta P [mmHg]
y(2) = Plv; % Left Ventricle pressure [mmHg]
y(3) = Psa; % Systemic Aorta pressure [mmHg]
y(4) = DA;  % Aortial valve status
y(5) = DM;  % Mitral valve status
end