% Perform Symbolic analysis necessary for obtaining scaled states and
% measurements from those states.
clear all; close all; clc;

syms qao qr Vlv lambda Rsvr Cao Cr Ra Rm Rc Lc s E;
syms qaos qrs Vlvs lambdas v0s Qvads;

%% Linear system definition
A = [-1/(Rsvr*Cao), 1/(Rsvr*Cr), 0, 1/Lc;
     1/(Rsvr*Cao), -(1/(Rsvr*Cr) + 1/(Rm*Cr)), E/Rm, 0;
     0, 1/(Rm*Cr), -E/Rm, -1/Lc;
     -1/Cao, 0, E, -(Ra+Rc)/Lc];
B = [0; -E/Rm; E/Rm; -E];
C = [0; 0; 0; Ra];

%% Scaling
% Sx = diag([300, 50, 200, 2]);
Sx = diag([300, 50, 150, 1.4])    % Updated to this after seeing vad = on, unhealthy case. This covers all those state variations as well.
Sv0 = [5];
Sqvad = [5*1000/60];

As = (inv(Sx))*A*Sx;
Bs = (inv(Sx))*B*Sv0;
Cs = (inv(Sx))*C*Sqvad;

xs = [qaos; qrs; Vlvs; lambdas];

fprintf('The scaled state equations are \n');
As*xs + Bs*v0s + Cs*Qvads
pretty(As*xs + Bs*v0s + Cs*Qvads);

%% Measurements in terms of scaled states (actual measurements)
fA = [0,0,0,1/Lc]*Sx*xs;
Plv = [0,0,E,0]*Sx*xs - E*Sv0*v0s;
Psa = [0,0,E,-Ra/Lc]*Sx*xs - E*Sv0*v0s + Ra*Sqvad*Qvads;
Qvad = Sqvad*Qvads;

fprintf('1. Delta pressure \n');
pretty(Plv - Psa);

fprintf('2. Plv pressure in terms of scaled states \n');
pretty(Plv);

fprintf('3. Psa in terms of scaled states \n ');
pretty(Psa);

fprintf('4. fA in terms of scaled states \n');
pretty(fA);

fprintf('5. Qvad in terms of scaled states \n');
pretty(Qvad);

fprintf('6. qr/cr in terms of scaled states \n');
pretty([0,1/Cr,0,0]*Sx*xs);

%% Results (Diary)
% The scaled state equations are 
% /                   lambdas     qaos        qrs                      \
% |                   ------- - -------- + ---------                   |
% |                    150 Lc   Cao Rsvr   6 Cr Rsvr                   |
% |                                                                    |
% |         6 qaos        /   1        1    \   4 E Vlvs   E v0s       |
% |        -------- - qrs | ----- + ------- | + -------- - -----       |
% |        Cao Rsvr       \ Cr Rm   Cr Rsvr /      Rm      10 Rm       |
% |                                                                    |
% |                   qrs     lambdas   E Vlvs   E v0s                 |
% |                 ------- - ------- - ------ + -----                 |
% |                 4 Cr Rm    100 Lc     Rm     40 Rm                 |
% |                                                                    |
% | 125 Qvads Ra   5 E v0s   150 qaos                lambdas (Ra + Rc) |
% | ------------ - ------- - -------- + 100 E Vlvs - ----------------- |
% \       3           2         Cao                          Lc        /
% 
% 1. Delta pressure 
% 2 Ra lambdas   250 Qvads Ra
% ------------ - ------------
%      Lc              3
% 
% 2. Plv pressure in terms of scaled states 
% 200 E Vlvs - 5 E v0s
% 
% 3. Psa in terms of scaled states 
%  250 Qvads Ra                          2 Ra lambdas
% ------------ - 5 E v0s + 200 E Vlvs - ------------
%       3                                    Lc
% 
% 4. fA in terms of scaled states 
% 2 lambdas
% ---------
%     Lc
% 
% 5. Qvad in terms of scaled states 
% 250 Qvads
% ---------
%     3
% 
% 6. qr/cr in terms of scaled states 
% 50 qrs
% ------
%   Cr