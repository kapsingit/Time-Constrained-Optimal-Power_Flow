function [y] = hess_func(xt,lamt,mut,T,d2ft,Sflist,Stlist,dSf_dValist,dSf_dVmlist,dSt_dVmlist,dSt_dValist)

[bus_data,gen_data,branch_data, gen_cost,baseMVA] = network_info;

[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C,ANGMIN, ANGMAX ] = idfor_branch;

[Ybus, Yf, Yt] = make_ybus;

il = find(branch_data(:,RATE_A)~=0);
nl2 = length(il);
nb = size(bus_data, 1);

y = [];

f2= branch_data(il, F_BUS);    %% list of "from" buses
t2 = branch_data(il, T_BUS);    %% list of "to" buses
Cf = sparse(1:nl2, f2, ones(nl2, 1), nl2, nb);     %% connection matrix for line & from buses
Ct = sparse(1:nl2, t2, ones(nl2, 1), nl2, nb);     %% connection matrix for line & to buses

for i= 1:T
    
    x = xt([1:20]+(i-1)*20);
    lam = lamt([1:11]+(i-1)*11);
    mu = mut([1:34]+(i-1)*34);
    d2f = d2ft(([1:20]+(i-1)*20):([1:20]+(i-1)*20));
    
    Va = x(1:5,:);
    Vm = x(6:10,:);
    Pg = x(11:15,:);
    Qg = x(16:20,:);

    V = Vm .* exp(1j * Va);
    lamP = lam(1:5);
    lamQ = lam(6:10);
    
    nxtra = length(x) - 2*nb;
    
    [Gpaa, Gpav, Gpva, Gpvv] = d2Sbus_dV2(Ybus, V, lamP);
    [Gqaa, Gqav, Gqva, Gqvv] = d2Sbus_dV2(Ybus, V, lamQ);
    
    d2G = [
    real([Gpaa Gpav; Gpva Gpvv]) + imag([Gqaa Gqav; Gqva Gqvv]) sparse(2*nb, nxtra);
    sparse(nxtra, 2*nb + nxtra)
    ];
    
    muF = mu(1:2);
    muT = mu(3:4);
    
    [Hfaa, Hfav, Hfva, Hfvv] = d2ASbr_dV2(dSf_dValist(i).value, dSf_dVmlist(i).value, Sflist(i).value, Cf, Yf(il,:), V, muF);
    [Htaa, Htav, Htva, Htvv] = d2ASbr_dV2(dSt_dValist(i).value, dSt_dVmlist(i).value, Stlist(i).value, Ct, Yt(il,:), V, muT);

    d2H = [
    [Hfaa Hfav; Hfva Hfvv] + [Htaa Htav; Htva Htvv] sparse(2*nb, nxtra);
    sparse(nxtra, 2*nb + nxtra)
    ];
    
    y = blkdiag(y,d2f + d2G + d2H);
end




