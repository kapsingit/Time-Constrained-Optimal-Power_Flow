function [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C,ANGMIN, ANGMAX ] = idfor_branch

F_BUS       = 1;    %% f, from bus number
T_BUS       = 2;    %% t, to bus number
BR_R        = 3;    %% r, resistance (p.u.)
BR_X        = 4;    %% x, reactance (p.u.)
BR_B        = 5;    %% b, total line charging susceptance (p.u.)
RATE_A      = 6;    %% rateA, MVA rating A (long term rating)
RATE_B      = 7;    %% rateB, MVA rating B (short term rating)
RATE_C      = 8;    %% rateC, MVA rating C (emergency rating)
ANGMIN      = 9;   %% minimum angle difference, angle(Vf) - angle(Vt) (degrees)
ANGMAX      = 10;   %% maximum angle difference, angle(Vf) - angle(Vt) (degrees)