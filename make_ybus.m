function [Ybus, Yf, Yt] = make_ybus

[bus_data,gen_data,branch_data, gen_cost,baseMVA] = network_info;

[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, VM, ...
    VA, BASE_KV, VMAX, VMIN] = idfor_bus;

[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C,ANGMIN, ANGMAX ] = idfor_branch;


Ys = 1 ./ (branch_data(:, BR_R) + 1j * branch_data(:, BR_X));  %% series admittance, 1 is status in original matpower

Bc = 1 .* branch_data(:, BR_B);

nb = size(bus_data, 1);          %% number of buses
nl = size(branch_data, 1);       %% number of lines

Ytt = Ys + 1j*Bc/2;
Yff = Ytt;
Yft = - Ys;
Ytf = - Ys;

Ysh = (bus_data(:, GS) + 1j * bus_data(:, BS)) / baseMVA;

f = branch_data(:, F_BUS);                           %% list of "from" buses
t = branch_data(:, T_BUS);                           %% list of "to" buses
Cf = sparse(1:nl, f, ones(nl, 1), nl, nb);      %% connection matrix for line & from buses
Ct = sparse(1:nl, t, ones(nl, 1), nl, nb);      %% connection matrix for line & to buses

i = [1:nl; 1:nl]'; 
Yf = sparse(i, [f; t], [Yff; Yft], nl, nb);
Yt = sparse(i, [f; t], [Ytf; Ytt], nl, nb);
Ybus = Cf' * Yf + Ct' * Yt + ...                %% branch admittances
        sparse(1:nb, 1:nb, Ysh, nb, nb); 