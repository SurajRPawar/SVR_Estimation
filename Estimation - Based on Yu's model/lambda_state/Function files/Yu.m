% Function file to simulate scaled Yu model. 
% Inputs - Current states, Current Elastance, Current QVAD, End diastolic
% volume
function [xdot, y] = Yu(x,E,Qvad,V0)
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
fA = lambda/Lc;
Plv = E*(Vlv - V0);
Psa = Qvad*Ra + Plv - Ra*fA;
Pla = qr/Cr; 
Pst = qao/Cao;

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
qaodot = fA - (1/Rsvr)*(Pst - Pla);
qrdot  = (1/Rsvr)*(Pst - Pla) - (1/Rm)*(Pla - Plv);
Vlvdot = (1/Rm)*(Pla - Plv) - fA;
lambdadot = Plv - Ra*(fA - Qvad) - Rc*fA - Pst;

xdot = [qaodot; qrdot; Vlvdot; lambdadot];

% Outputs
y(1) = Plv - Psa; % Delta P [mmHg]
y(2) = Plv; % Left Ventricle pressure [mmHg]
y(3) = Psa; % Systemic Aorta pressure [mmHg]
y(4) = DA;  % Aortial valve status
y(5) = DM;  % Mitral valve status
y(6) = fA;
end

%% Symbolic analysis
% syms qao qr Vlv Vo lambda Rsvr Cao Cr Lc Ra Rm E w Qvad Rc;
% x = [qao; qr; Vlv; lambda; Rsvr; Cao];
% u = [Vo];
% f = [lambda/Lc - (1/Rsvr)*(qao/Cao - qr/Cr) + w;
%     (1/Rsvr)*(qao/Cao - qr/Cr) - (1/Rm)*(qr/Cr - E*(Vlv-Vo)) + w;
%     (1/Rm)*(qr/Cr - E*(Vlv-Vo)) - Qvad - 1*(lambda/Lc - Qvad) + w;
%     1*(E*(Vlv-Vo) - 1*Ra*(lambda/Lc-Qvad)) - Rc*lambda/Lc - 1*qao/Cao + w;
%     0;
%     0;
%     0];
% g = [Ra*(lambda/Lc - Qvad)];
% 
% Fsym = jacobian(f,x)
% Hsym = jacobian(g,x)
% Gamma = jacobian(f,w)
