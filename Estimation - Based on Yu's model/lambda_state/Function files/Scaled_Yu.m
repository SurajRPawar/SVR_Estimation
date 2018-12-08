% Function file to simulate scaled Yu model. 
% Inputs - Current states, Current Elastance, Current QVAD, End diastolic
% volume
function [xdot, y] = Scaled_Yu(x,E,Qvad,V0)
% x = States [qao; qr; Vlv; lambda]
% E = Emax*e_t_lv(i), elastance
% Qvad = VAD flow [mL/sec]
% V0 = End diastolic volume in LV [mL]

global Rc Lc Cao Rsvr Cr Ra Rm amp Rv;

% Define states
qao     = x(1); 
qr      = x(2);
Vlv     = x(3);
lambda  = x(4);

% Calculate required variables
fA = 7*lambda/(5*Lc);
Plv = E*(150*Vlv - 5*V0);
Psa = 250*Qvad*Ra/3 + Plv - Ra*fA;
Pla = 50*qr/Cr; 

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

% State equations
% qaodot    = lambda/(150*Lc) - (1/Rsvr)*(qao/Cao - qr/(6*Cr));
% qrdot     = 6*qao/(Rsvr*Cao) - qr*(1/(Cr*Rm) + 1/(Cr*Rsvr)) + 4*E*Vlv/Rm - E*V0/(10*Rm);
% Vlvdot    = qr/(4*Cr*Rm) - lambda/(100*Lc) - E*Vlv/Rm + E*V0/(40*Rm);
% lambdadot = 125*Qvad*Ra/3 - 5*E*V0/2 - 150*qao/Cao + 100*E*Vlv - lambda*(Ra+Rc)/Lc;

% Below are new equations, updated based on new scaling matrix
qaodot    = 7*lambda/(1500*Lc) - qao/(Cao*Rsvr) + qr/(6*Cr*Rsvr);
qrdot     = (6*qao)/(Cao*Rsvr) - qr*(1/(Cr*Rm) + 1/(Cr*Rsvr)) + (3*E*Vlv)/Rm - (E*V0)/(10*Rm);
Vlvdot    = qr/(3*Cr*Rm) - 7*lambda/(750*Lc) - (E*Vlv)/Rm + (E*V0)/(30*Rm);
lambdadot = (1250*Qvad*Ra)/21 - 25*E*V0/7 - (1500*qao)/(7*Cao) + 750*E*Vlv/7 - (lambda*(Ra + Rc))/Lc;
     
xdot = [qaodot; qrdot; Vlvdot; lambdadot];

% Outputs
y(1) = Plv - Psa; % Delta P [mmHg]
y(2) = Plv; % Left Ventricle pressure [mmHg]
y(3) = Psa; % Systemic Aorta pressure [mmHg]
y(4) = DA;  % Aortial valve status
y(5) = DM;  % Mitral valve status
end

