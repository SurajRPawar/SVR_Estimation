% Version 5 - Estimating only Rsvr
function [fx_est, Pdot, y] = Yu_est_v1(xhat,E,Qvad,V0,P,Q,DA,DM);
global Rc Lc Cao Rsvr Cr Ra Rm amp Rv;

% States and estimates
qaohat    = xhat(1);
qrhat     = xhat(2);
Vlvhat    = xhat(3);
fAhat     = xhat(4);
Rsvrhat   = xhat(5);

% Variables to calculate
Plv = E*(Vlvhat - V0);
Psa = Qvad*Ra + Plv - Ra*fAhat;
Pla = qrhat/Cr; 

% Toggle valve status

if DA == 1;
% if Plv > Psa
    % Ventricle Ejection stage, Aortic valve open (1)
%     DA = 1;
    Ra = Rv;
else
    % Ventricle filling stage, Aortic valve closed (0)
%     DA = 0;
    Ra = amp*Rv;
end

if DM == 1;
% if Pla > Plv
    % Atria ejection stage, Mitral valve open (1)
%     DM = 1;
    Rm = Rv;
else
%     DM = 0;
    Rm = amp*Rv;
end

% fx_est = [lambdahat/(150*Lc) - (1/Rsvrhat)*(qaohat/Cao - qrhat/(6*Cr));
%           6*qaohat/(Rsvrhat*Cao) - qrhat*(1/(Cr*Rm) + 1/(Cr*Rsvrhat)) + 4*E*Vlvhat/Rm - E*V0s/(10*Rm);
%           qrhat/(4*Cr*Rm) - lambdahat/(100*Lc) - E*Vlvhat/Rm + E*V0s/(40*Rm);
%           125*Qvad*Ra/3 - 5*E*V0s/2 - 150*qaohat/Cao + 100*E*Vlvhat - lambdahat*(Ra+Rc)/Lc;
%           0];

% New equations with new scaling matrix
qaodot    = fAhat - qaohat/(Cao*Rsvrhat) + qrhat/(Cr*Rsvrhat);
qrdot     = (1/Rsvrhat)*(qaohat/Cao - qrhat/Cr) - (1/Rm)*(qrhat/Cr - Plv);
Vlvdot    = (1/Rm)*(qrhat/Cr - Plv) - fAhat;
fAdot     = (1/Lc)*(Plv - Ra*(fAhat - Qvad) - Rc*fAhat - qaohat/Cao);
fx_est    = [qaodot; qrdot; Vlvdot; fAdot; 0];

% New F matrix, using new scaling matrix Sx
F = ...
[ -1/(Cao*Rsvrhat),               1/(Cr*Rsvrhat),     0,             1, qaohat/(Cao*Rsvrhat^2) - qrhat/(Cr*Rsvrhat^2);
   1/(Cao*Rsvrhat), - 1/(Cr*Rm) - 1/(Cr*Rsvrhat),  E/Rm,             0,            -(qaohat/Cao - qrhat/Cr)/Rsvrhat^2;
                 0,                    1/(Cr*Rm), -E/Rm,            -1,                                             0;
       -1/(Cao*Lc),                            0,  E/Lc, -(Ra + Rc)/Lc,                                             0;
                 0,                            0,     0,             0,                                             0];
             
Pdot = F*P + P*F.' + Q;

% Outputs
y(1) = Plv - Psa; % Delta P [mmHg]
y(2) = Plv; % Left Ventricle pressure [mmHg]
y(3) = Psa; % Systemic Aorta pressure [mmHg]
% y(4) = DA;  % Aortial valve status
% y(5) = DM;  % Mitral valve status

%% Symbolic analysis
% syms qaohat qrhat Vlvhat V0 fAhat Rsvrhat Cao Cr Lc Ra Rm E w Qvad Rc;
% syms qaodot qrdot Vlvdot fAdot;
% 
% x = [qaohat; qrhat; Vlvhat; fAhat; Rsvrhat];
% u = [V0];
% 
% Plv = E*(Vlvhat - V0);
% qaodot    = fAhat - qaohat/(Cao*Rsvrhat) + qrhat/(Cr*Rsvrhat) + w;
% qrdot     = (1/Rsvrhat)*(qaohat/Cao - qrhat/Cr) - (1/Rm)*(qrhat/Cr - Plv) + w;
% Vlvdot    = (1/Rm)*(qrhat/Cr - Plv) - fAhat + w;
% fAdot     = (1/Lc)*(Plv - Ra*(fAhat - Qvad) - Rc*fAhat - qaohat/Cao) + w;
% f    = [qaodot; qrdot; Vlvdot; fAdot; 0];
% 
% g = [Ra*(fAhat - Qvad)];
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

