function [x, xmin, xmax] = decesion_variables(T)
% setting up total number of problem variables;
% each bus will have two variables voltage magnitude, voltage angle,
% each generator will have active power injection and reactive power injection
% variables are arranged like [Va,Vm,Pg,Qg] each with five entries so in
% total 20 variables for 5 bus system at one time

% For 24 hour long TCOPF in 1 hour time step 
% number of variables = 24 x 20 = 480
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, VM, ...
    VA, BASE_KV, VMAX, VMIN] = idfor_bus;

[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, RAMP] = idfor_gen;

[bus_data,gen_data,branch_data, gen_cost,baseMVA] = network_info;

nb = size(bus_data, 1); 

Va_max = Inf(nb,1);
Va_min = -Va_max;
refs  = bus_data(:,BUS_TYPE)==3;
Va_max(refs) = bus_data(refs,VA); % Minimum limit for voltage angle
Va_min(refs) = bus_data(refs,VA); % Maximum limit for voltage angle
Va = bus_data(:,VA); % Initial voltage angle
Vm_max = 1.1*ones(nb,1); % Voltage magnitude maximum limit
Vm_min = 0.9*ones(nb,1); % Voltage magnitude minimum limit
Vm = bus_data(:,VM); % Voltage magnitude of bus
Pgmax = gen_data(:,PMAX)/100;
Pgmin = gen_data(:,PMIN)/100;
Pg = (Pgmax+Pgmin)/2;
Qgmax = gen_data(:,QMAX)/100;
Qgmin = gen_data(:,QMIN)/100;
Qg = (Qgmax + Qgmin)/2;
x = [Va;Vm;Pg;Qg]; % vector of the decision variables
xmin = [Va_min;Vm_min;Pgmin;Qgmin];
xmax = [Va_max;Vm_max;Pgmax;Qgmax];

%% For T time steps we put the same initial values
x = repmat(x,T,1);
xmax = repmat(xmax, T, 1);
xmin = repmat(xmin, T, 1);

