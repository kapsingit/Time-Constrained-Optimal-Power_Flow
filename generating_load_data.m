function [load_data_p, load_data_q] = generating_load_data(T)


[bus_data,gen_data,branch_data, gen_cost,baseMVA] = network_info;

[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, VM, ...
    VA, BASE_KV, VMAX, VMIN] = idfor_bus;

%a = 0.5 + 0.5*rand(1,T); % Varying load radomly between 50% to 100% to generate synthetic load pattern
a =[ 1, 0.5 + 0.5*rand(1,T-1)];
load_data_p = zeros(5,T);
load_data_q = zeros(5,T);

for i = 1:T
    load_data_p(:,i) = bus_data(:,PD)*a(i);
    load_data_q(:,i) = bus_data(:,QD)*a(i);
end