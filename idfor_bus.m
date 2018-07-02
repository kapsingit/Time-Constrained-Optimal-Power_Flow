function [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, VM, ...
    VA, BASE_KV, VMAX, VMIN] = idfor_bus

%% define bus types
PQ      = 1;
PV      = 2;
REF     = 3;
NONE    = 4;

%% define the indices
BUS_I       = 1;    %% bus number (1 to 29997)
BUS_TYPE    = 2;    %% bus type (1 - PQ bus, 2 - PV bus, 3 - reference bus, 4 - isolated bus)
PD          = 3;    %% Pd, real power demand (MW)
QD          = 4;    %% Qd, reactive power demand (MVAr)
GS          = 5;    %% Gs, shunt conductance (MW at V = 1.0 p.u.)
BS          = 6;    %% Bs, shunt susceptance (MVAr at V = 1.0 p.u.)
VM          = 7;    %% Vm, voltage magnitude (p.u.)
VA          = 8;    %% Va, voltage angle (degrees)
BASE_KV     = 9;   %% baseKV, base voltage (kV)
VMAX        = 10;   %% maxVm, maximum voltage magnitude (p.u.)      (not in PTI format)
VMIN        = 11;   %% minVm, minimum voltage magnitude (p.u.)      (not in PTI format)

