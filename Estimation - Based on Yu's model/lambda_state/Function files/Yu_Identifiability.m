% Identifiability analysis for Yu- inspired model

close all; clear all; clc;

syms s qao qr Vlv lambda E V0 Rsvr Cao Rc Lc Cr Ra Rm Qvad;

x = [qao; qr; Vlv; lambda];
u1 = [Qvad];
u2 = [V0];
f = [lambda/Lc - (1/Rsvr)*(qao/Cao - qr/Cr);
     (1/Rsvr)*(qao/Cao - qr/Cr) - (1/Rm)*(qr/Cr - E*(Vlv-V0));
     (1/Rm)*(qr/Cr - E*(Vlv - V0)) - lambda/Lc;
     E*(Vlv - V0) - Ra*(lambda/Lc-Qvad) - Rc*lambda/Lc - qao/Cao];
g = [Ra*(lambda/Lc - Qvad)];

A = jacobian(f,x);
B = jacobian(f,u1);
C = jacobian(g,x);
D = jacobian(g,u1);

phi = inv(s*eye(length(x)) - A);
TF = C*phi*B + D;
[num, den] = numden(TF);
collect(num,s)
collect(den,s)
% pretty(TF);