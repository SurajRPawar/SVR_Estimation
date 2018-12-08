close all; clear all; clc;

% syms qao qr Vlv lambda Rsvr Cao Cr Qvad Ra Rm Rc Lc;
syms qaohat qrhat Vlvhat lambdahat Rsvrhat Caohat Crhat Qvad Ra Rm Rc Lc x5hat x6hat x7hat E V0s;

x = [qaohat; qrhat; Vlvhat; lambdahat; x5hat; x6hat; x7hat];
% x5 = 1/(Rsvr*Cr)
% x6 = 1/(Rsvr*Cao)
% x7 = 1/Cao

% For 3 parameters estimation
% f = [lambdahat/(150*Lc) - qaohat*x6hat + qrhat*x5hat/6; 
%      6*qaohat*x6hat - qrhat/(Crhat*Rm) - qrhat*x5hat + 4*E*Vlvhat/Rm - E*V0/(10*Rm);
%      qrhat/(4*Crhat*Rm) - lambdahat/(100*Lc) - E*Vlvhat/Rm + E*V0/(40*Rm);
%      125*Qvad*Ra/3 - 5*E*V0/2 - 150*qaohat*x7hat + 100*E*Vlvhat - lambdahat*(Ra+Rc)/Lc;
%      0;
%      0;
%      0];

% For 2 parameter estimation with lumped parameters
% f = [lambdahat/(150*Lc) - qaohat*x5hat*x6hat*Cr + qrhat*x5hat/6; 
%      6*qaohat*x5hat*x6hat*Cr - qrhat/(Cr*Rm) - qrhat*x5hat + 4*E*Vlvhat/Rm - E*V0s/(10*Rm);
%      qrhat/(4*Cr*Rm) - lambdahat/(100*Lc) - E*Vlvhat/Rm + E*V0s/(40*Rm);
%      125*Qvad*Ra/3 - 5*E*V0s/2 - 150*qaohat*x6hat + 100*E*Vlvhat - lambdahat*(Ra+Rc)/Lc;
%      0;
%      0];

% For 3 parameters independent.

f = [lambdahat/(150*Lc) - qaohat*x5hat*x6hat + qrhat*x5hat*x7hat/6; 
          6*qaohat*x5hat*x6hat - qrhat*(x7hat/Rm + x5hat*x7hat) + 4*E*Vlvhat/Rm - E*V0s/(10*Rm);
          qrhat*x7hat/(4*Rm) - lambdahat/(100*Lc) - E*Vlvhat/Rm + E*V0s/(40*Rm);
          125*Qvad*Ra/3 - 5*E*V0s/2 - 150*x6hat*qaohat + 100*E*Vlvhat - lambdahat*(Ra+Rc)/Lc;
          0;
          0;
          0];
F = jacobian(f,x)
% pretty(F);

% Without augmenting, this is F
% /       1             1                         1    \
% | - --------,     ---------,         0,      ------  |
% |   Cao Rsvr      6 Cr Rsvr                  150 Lc  |
% |                                                    |
% |      6          1        1         4E               |
% |  --------,  - ----- - -------,   -----,      0     |
% |  Cao Rsvr     Cr Rm   Cr Rsvr       Rm             |
% |                                                    |
% |                     1               E          1   |
% |      0,          -------,      - ------,  - ------ |
% |                  4 Cr Rm            Rm      100 Lc |
% |                                                    |
% |      150                                   Ra + Rc |
% |    - ---,           0,          100E,    - ------- |
% \      Cao                                      Lc   /

% After state augmentation

% F =
%  
% [     -x6hat,                x5hat/6,        0,    1/(150*Lc), qrhat/6,  -qaohat,           0]
% [    6*x6hat, - x5hat - 1/(Crhat*Rm), (4*E)/Rm,             0,  -qrhat, 6*qaohat,           0]
% [          0,         1/(4*Crhat*Rm),    -E/Rm,   -1/(100*Lc),       0,        0,           0]
% [ -150*x7hat,                      0,    100*E, -(Ra + Rc)/Lc,       0,        0, -150*qaohat]
% [          0,                      0,        0,             0,       0,        0,           0]
% [          0,                      0,        0,             0,       0,        0,           0]
% [          0,                      0,        0,             0,       0,        0,           0]


% Pretty format - 
% /                    x5hat                    1      qrhat                        \
% |   -x6hat,          -----,         0,     ------,   -----,  -qaohat,      0      |
% |                      6                   150 Lc      6                          |
% |                                                                                 |
% |                           1      4 E                                            |
% |   6 x6hat,  - x5hat - --------,  ---,      0,     -qrhat, 6 qaohat,      0      |
% |                       Crhat Rm    Rm                                            |
% |                                                                                 |
% |                      1              E        1                                  |
% |      0,         ----------,      - --,  - ------,    0,       0,         0      |
% |                 4 Crhat Rm         Rm     100 Lc                                |
% |                                                                                 |
% |                                          Ra + Rc                                |
% | -150 x7hat,          0,         100 E, - -------,    0,       0,    -150 qaohat |
% |                                             Lc                                  |
% |                                                                                 |
% |      0,              0,           0,       0,        0,       0,         0      |
% |                                                                                 |
% |      0,              0,           0,       0,        0,       0,         0      |
% |                                                                                 |
% \      0,              0,           0,       0,        0,       0,         0      /

% For 2 parameter estimation

% F =
%  
% [  -Cr*x5hat*x6hat,             x5hat/6,        0,    1/(150*Lc), qrhat/6 - Cr*qaohat*x6hat,  -Cr*qaohat*x5hat]
% [ 6*Cr*x5hat*x6hat, - x5hat - 1/(Cr*Rm), (4*E)/Rm,             0, 6*Cr*qaohat*x6hat - qrhat, 6*Cr*qaohat*x5hat]
% [                0,         1/(4*Cr*Rm),    -E/Rm,   -1/(100*Lc),                         0,                 0]
% [       -150*x6hat,                   0,    100*E, -(Ra + Rc)/Lc,                         0,       -150*qaohat]
% [                0,                   0,        0,             0,                         0,                 0]
% [                0,                   0,        0,             0,                         0,                 0]
% 
% pretty(F)
% /                        x5hat                   1      qrhat                                       \
% |  -Cr x5hat x6hat,      -----,        0,     ------,   ----- - Cr qaohat x6hat,   -Cr qaohat x5hat |
% |                          6                  150 Lc      6                                         |
% |                                                                                                   |
% |                               1     4 E                                                           |
% | 6 Cr x5hat x6hat, - x5hat - -----,  ---,      0,     6 Cr qaohat x6hat - qrhat, 6 Cr qaohat x5hat |
% |                             Cr Rm    Rm                                                           |
% |                                                                                                   |
% |                          1             E        1                                                 |
% |         0,            -------,      - --,  - ------,             0,                     0         |
% |                       4 Cr Rm         Rm     100 Lc                                               |
% |                                                                                                   |
% |                                             Ra + Rc                                               |
% |    -150 x6hat,           0,        100 E, - -------,             0,                -150 qaohat    |
% |                                                Lc                                                 |
% |                                                                                                   |
% |         0,               0,          0,       0,                 0,                     0         |
% |                                                                                                   |
% \         0,               0,          0,       0,                 0,                     0         /