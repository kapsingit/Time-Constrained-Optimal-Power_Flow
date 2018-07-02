function [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, RAMP] = idfor_gen

%% define the indices
GEN_BUS     = 1;    %% bus number
PG          = 2;    %% Pg, real power output (MW)
QG          = 3;    %% Qg, reactive power output (MVAr)
QMAX        = 4;    %% Qmax, maximum reactive power output at Pmin (MVAr)
QMIN        = 5;    %% Qmin, minimum reactive power output at Pmin (MVAr)
VG          = 6;    %% Vg, voltage magnitude setpoint (p.u.)
MBASE       = 7;    %% mBase, total MVA base of this machine, defaults to baseMVA
GEN_STATUS  = 8;    %% status, 1 - machine in service, 0 - machine out of service
PMAX        = 9;    %% Pmax, maximum real power output (MW)
PMIN        = 10;   %% Pmin, minimum real power output (MW)
RAMP        = 11;
