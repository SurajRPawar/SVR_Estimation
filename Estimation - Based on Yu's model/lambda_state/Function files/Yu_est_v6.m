% Version 6 - Estimating only Rsvr. Same as version 5, but trying to use a
% switch like structure.

function [fx_est, Pdot, y] = Yu_est_v1(xhat,E,Qvad,V0s,P,Q);
global Rc Lc Cao Rsvr Cr Ra Rm amp Rv;

% States and estimates
qaohat    = xhat(1);
qrhat     = xhat(2);
Vlvhat    = xhat(3);
lambdahat = xhat(4);
Rsvrhat   = xhat(5);

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

fx_est = [fA/300 - (1/Rsvrhat)*(qaohat/Cao - qrhat/(6*Cr));
          6*qaohat/(Rsvrhat*Cao) - qrhat*(1/(Cr*Rm) + 1/(Cr*Rsvrhat)) + 4*E*Vlvhat/Rm - E*V0s/(10*Rm);
          qrhat/(4*Cr*Rm) - lambdahat/(100*Lc) - E*Vlvhat/Rm + E*V0s/(40*Rm);
          125*Qvad*Ra/3 - 5*E*V0s/2 - 150*qaohat/Cao + 100*E*Vlvhat - lambdahat*(Ra+Rc)/Lc;
          0];

F = ...
[ -1/(Cao*Rsvrhat),             1/(6*Cr*Rsvrhat),        0,    1/(150*Lc),             (qaohat/Cao - qrhat/(6*Cr))/Rsvrhat^2;
   6/(Cao*Rsvrhat), - 1/(Cr*Rm) - 1/(Cr*Rsvrhat), (4*E)/Rm,             0, qrhat/(Cr*Rsvrhat^2) - (6*qaohat)/(Cao*Rsvrhat^2);
                 0,                  1/(4*Cr*Rm),    -E/Rm,   -1/(100*Lc),                                                 0;
          -150/Cao,                            0,    100*E, -(Ra + Rc)/Lc,                                                 0;
                 0,                            0,        0,             0,                                                 0];
 
      
Pdot = F*P + P*F.' + Q;

% Outputs
y(1) = Plv - Psa; % Delta P [mmHg]
y(2) = Plv; % Left Ventricle pressure [mmHg]
y(3) = Psa; % Systemic Aorta pressure [mmHg]
y(4) = DA;  % Aortial valve status
y(5) = DM;  % Mitral valve status

%% Symbolic analysis
% syms qaohat qrhat Vlvhat V0s lambdahat Rsvrhat Cao Cr Lc Ra Rm E w Qvad Rc;
% x = [qaohat; qrhat; Vlvhat; lambdahat; Rsvrhat];
% u = [V0s];
% 
% f = ...
% [lambdahat/(150*Lc) - (1/Rsvrhat)*(qaohat/Cao - qrhat/(6*Cr));
% 6*qaohat/(Rsvrhat*Cao) - qrhat*(1/(Cr*Rm) + 1/(Cr*Rsvrhat)) + 4*E*Vlvhat/Rm - E*V0s/(10*Rm);
% qrhat/(4*Cr*Rm) - lambdahat/(100*Lc) - E*Vlvhat/Rm + E*V0s/(40*Rm);
% 125*Qvad*Ra/3 - 5*E*V0s/2 - 150*qaohat/Cao + 100*E*Vlvhat - lambdahat*(Ra+Rc)/Lc;
% 0];
% 
% g = [2*Ra*lambdahat/Lc - 250*Qvad*Ra/3];
% 
% Fsym = jacobian(f,x)
% Hsym = jacobian(g,x)
% Gamma = jacobian(f,w)

% Try using these gains, they gave best results with this version
% Qk = 0.6; % Process noise, Higher Q -> Higher transients, less confidence on model
% Pk = 0.3; % Initial error covariance. Higher P -> Lesser confidence in IC
% 
% Q = diag(Qk.*[1,1,1,1]);
% P = diag(Pk.*[1,1,1,1]);
% R = (0.5)^2;
end

