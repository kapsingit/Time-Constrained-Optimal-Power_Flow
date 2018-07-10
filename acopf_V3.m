clear;clc;

timesteps = 24;
[bus_data,gen_data,branch_data, gen_cost,baseMVA] = network_info;

load loaddata;

%% Clearing workspace and command window
%%  Making admittance matrix 

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

%% 
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN] = idfor_gen;

% setting up total number of problem variables;
% each bus will have two variables voltage magnitude, voltage angle,
% each generator will have active power injection and reactive power injection
% variables are arranged like [Va,Vm,Pg,Qg] each with five entries so in
% total 20 variables for 5 bus system

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

xt = repmat(x,timesteps,1);
xmint = repmat(xmin,timesteps,1);
xmaxt = repmat(xmax,timesteps,1);

ng = size(gen_data,1);
n_var = length(xt);
il = [1,6]; % lines 1 and 6 has the limits for line flow;
nl2 = length(il); % number of constrained lines

%%
ht = [];gt = []; dht = []; dgt = []; gnt = []; hnt = [];
for i = 1:timesteps
    x = xt([1:20]+(i-1)*20);
    [h, g, dh, dg,gn,hn,dSf_dVa,dSf_dVm,dSt_dVm,dSt_dVa,Sf,St] = gh_fcn1(load_data_p(:,i),load_data_q(:,i),x,Ybus,bus_data,gen_data,branch_data,il,Yf,Yt,baseMVA,xmax,xmin);
    ht = [ht;h]; gt = [gt;g]; dht = blkdiag(dht,dh); dgt = blkdiag(dgt,dg); gnt = [gnt;gn]; hnt = [hnt;hn];
    Sflist(i).value = Sf;
    Stlist(i).value = St;
    dSf_dValist(i).value = dSf_dVa;
    dSf_dVmlist(i).value = dSf_dVm;
    dSt_dVmlist(i).value = dSt_dVm;
    dSt_dValist(i).value = dSt_dVa;
end
neq = size(gt, 1);           %% number of equality constraints
niq = size(ht, 1);           %% number of inequality constraints
neqnln = size(gnt, 1);       %% number of nonlinear equality constraints
niqnln = size(hnt, 1);       %% number of nonlinear inequality constraints

%% Now finding the objective function, its derivative and hessian 

d2f = sparse(1:20,1:20,0); % second derivative of the objective function

%% initialize gamma, lam, mu, z, e
z0 = 1;
gamma = 1;                  
lamt = zeros(neq, 1);
z   = z0 * ones(niq, 1);
mut  = z;
k = find(ht < -z0);
z(k) = -ht(k);
e = ones(niq, 1);

%% 

f2= branch_data(il, F_BUS);    %% list of "from" buses
t2 = branch_data(il, T_BUS);    %% list of "to" buses
Cf = sparse(1:nl2, f2, ones(nl2, 1), nl2, nb);     %% connection matrix for line & from buses
Ct = sparse(1:nl2, t2, ones(nl2, 1), nl2, nb);     %% connection matrix for line & to buses

%  
funt = 0;dft = [];
for i = 1:timesteps
    x = xt([1:20]+(i-1)*20);
    [fun, df] = f_fcn1(x,gen_cost,baseMVA); 
    funt = funt+fun;
    dft = [dft;df];
end

%%
kk = 1;
cond_OPF = [];
while kk <= 50
    zinvdiag = sparse(1:niq, 1:niq, 1 ./ z, niq, niq);
    mudiag = sparse(1:niq, 1:niq, mut, niq, niq);
    dh_zinv = dht * zinvdiag;
    Lx = dft + dgt * lamt + dht * mut;
    Lxxt = [];
    for i = 1:timesteps
        x = xt([1:20]+(i-1)*20);
        Sf = Sflist(i).value;
        St = Stlist(i).value;
        dSf_dVa = dSf_dValist(i).value;
        dSf_dVm = dSf_dVmlist(i).value;
        dSt_dVa = dSt_dValist(i).value;
        dSt_dVm = dSt_dVmlist(i).value;
        lam = lamt(1:length(lamt)/timesteps+(i-1)*length(lamt)/timesteps);
        mu = mut(1:length(mut)/timesteps+(i-1)*length(mut)/timesteps);
        Lxx = hess_fcn1(x,lam,mu,nb, Ybus,Yf,Yt,Cf,Ct,Sf,St,d2f,dSf_dVa,dSf_dVm,dSt_dVm,dSt_dVa);
        Lxxt = blkdiag(Lxxt,Lxx);
    end
    
    M = Lxxt + dh_zinv * mudiag * dht';
    N = Lx + dh_zinv * (mudiag * ht + gamma * e);
    W = [M dgt;dgt' sparse(neq, neq)];
    B = [-N; -gt];

    dxdlam = W \ B;
    
    %dxd(i,:) = dxdlam';
    nx = size(xt,1);
    dx = dxdlam(1:nx, 1);
    dlam = dxdlam(nx+(1:neq), 1);
    dz = -ht - z - dht' * dx;
    dmu = -mut + zinvdiag *(gamma*e - mudiag * dz);
    
    xi = 0.99995;
    k = find(dz < 0);
    if isempty(k)
        alphap = 1;
    else
        alphap = min( [xi * min(z(k) ./ -dz(k)) 1] );
    end
    
    k = find(dmu < 0);
    if isempty(k)
        alphad = 1;
    else
        alphad = min( [xi * min(mut(k) ./ -dmu(k)) 1] );
    end
    
    xt = xt + alphap * dx;
    z = z + alphap * dz;
    lamt = lamt + alphad * dlam;
    mut  = mut  + alphad * dmu;
    if niq > 0
        gamma = 0.1* (z' * mut) / niq;
    end
    funt = 0;dft = [];
    for i = 1:timesteps
    x = xt([1:20]+(i-1)*20);
    [fun, df] = f_fcn1(x,gen_cost,baseMVA); 
    funt = funt+fun;
    dft = [dft;df];
    end                 %% cost
    ht = [];gt = []; dht = []; dgt = []; gnt = []; hnt = [];
    for i = 1:timesteps
    x = xt([1:20]+(i-1)*20);
    [h, g, dh, dg,gn,hn,dSf_dVa,dSf_dVm,dSt_dVm,dSt_dVa,Sf,St] = gh_fcn1(load_data_p(:,i),load_data_q(:,i),x,Ybus,bus_data,gen_data,branch_data,il,Yf,Yt,baseMVA,xmax,xmin);
    ht = [ht;h]; gt = [gt;g]; dht = blkdiag(dht,dh); dgt = blkdiag(dgt,dg); gnt = [gnt;gn]; hnt = [hnt;hn];
    Sflist(i).value = Sf;
    Stlist(i).value = St;
    dSf_dValist(i).value = dSf_dVa;
    dSf_dVmlist(i).value = dSf_dVm;
    dSt_dVmlist(i).value = dSt_dVm;
    dSt_dValist(i).value = dSt_dVa;
    end
    solu(kk,:) = xt;
    kk = kk+1;
end

%%
sol = reshape(xt,20,timesteps);
PG = sol(11:15,:);
QG = sol(11:15,:);
figure('color',[1,1,1]);
plot(sum(PG*100));
hold
plot(sum(load_data_p));
legend('Gen','load')
xlabel('Time (in sec)');
ylabel('Power (MW)');
grid on; box on;
figure('color',[1,1,1]);
plot(solu(2:end,:)-solu(1:end-1,:));
xlabel('Iteration');
ylabel('Error magnitude');
grid on; box on;
%%

function [h, g, dh, dg,gn,hn,dSf_dVa,dSf_dVm,dSt_dVm,dSt_dVa,Sf,St] = gh_fcn1(pload,qload,x,Ybus,bus_data,gen_data,branch_data,il,Yf,Yt,baseMVA,xmax,xmin)

[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN] = idfor_gen;

[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, VM, ...
    VA, BASE_KV, VMAX, VMIN] = idfor_bus;

ng = size(gen_data,1);
nb = size(bus_data,1);

Pg = x(11:15,:);
Qg = x(16:20,:);

gen_buses = gen_data(:,GEN_BUS);

Cg = sparse(gen_buses, (1:ng)', 1, nb, ng); % generator connection matrix
%Cg = [1 1 0 0 0;0 0 0 0 0;0 0 1 0 0;0 0 0 1 0;0 0 0 0 1];
Sbusg = Cg*(Pg + 1j*Qg);
Sload = (pload+1j*qload)/baseMVA;

Sbus = Sbusg - Sload;

Vm = x(6:10,:);
Va = x(1:5,:);

V = Vm .* exp(1j * Va);
nl2 = length(il);

[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C,ANGMIN, ANGMAX ] = idfor_branch;

mis = V .* conj(Ybus * V) - Sbus; % calculating mismatch

gn = [ real(mis);            %% active power mismatch for all buses
      imag(mis) ];          %% reactive power mismatch for all buses

flow_max = (branch_data(il, RATE_A)/baseMVA).^2;
Sf = V(branch_data(il, F_BUS)) .* conj(Yf(il,:) * V);  %% complex power injected at "from" bus (p.u.)
St = V(branch_data(il, T_BUS)) .* conj(Yt(il,:) * V);  %% complex power injected at "to" bus (p.u.)

hn = [ Sf .* conj(Sf) - flow_max;      %% branch apparent power limits (from bus)
      St .* conj(St) - flow_max ];  %% branch apparent power limits (to bus)

% Finding the derivatives of the constraints

n = length(V);
ng = size(gen_data,1);
n_var = length(x);
Ibus = Ybus * V;
diagV       = sparse(1:n, 1:n, V, n, n);
diagIbus    = sparse(1:n, 1:n, Ibus, n, n);
diagVnorm   = sparse(1:n, 1:n, V./abs(V), n, n);
dSbus_dVm = diagV * conj(Ybus * diagVnorm) + conj(diagIbus) * diagVnorm;
dSbus_dVa = 1j * diagV * conj(diagIbus - Ybus * diagV);
neg_Cg = sparse(gen_data(:, GEN_BUS), 1:ng, -1, nb, ng);
dgn = sparse(2*nb, n_var);
dgn(:, [1:5 6:10 11:15 16:20]) = [
    real([dSbus_dVa dSbus_dVm]) neg_Cg sparse(nb, ng);  %% P mismatch w.r.t Va, Vm, Pg, Qg
    imag([dSbus_dVa dSbus_dVm]) sparse(nb, ng) neg_Cg;  %% Q mismatch w.r.t Va, Vm, Pg, Qg
];
dgn = dgn';

f1 = [1,4];
t1 = [2,5];
nl1 = length(f1);
If = Yf(il,:)* V;
It = Yt(il,:)* V;
Vnorm = V ./ abs(V);
diagVf      = sparse(1:nl1, 1:nl1, V(f1), nl1, nl1);
diagIf      = sparse(1:nl1, 1:nl1, If, nl1, nl1);
diagVt      = sparse(1:nl1, 1:nl1, V(t1), nl1, nl1);
diagIt      = sparse(1:nl1, 1:nl1, It, nl1, nl1);
diagV       = sparse(1:nb, 1:nb, V, nb, nb);
diagVnorm   = sparse(1:nb, 1:nb, Vnorm, nb, nb);
%     
dSf_dVa = 1j * (conj(diagIf) * sparse(1:nl1, f1, V(f1), nl1, nb) - diagVf * conj(Yf(il,:) * diagV));
dSf_dVm = diagVf * conj(Yf(il,:) * diagVnorm) + conj(diagIf) * sparse(1:nl1, f1, Vnorm(f1), nl1, nb);
dSt_dVa = 1j * (conj(diagIt) * sparse(1:nl1, t1, V(t1), nl1, nb) - diagVt * conj(Yt(il,:) * diagV));
dSt_dVm = diagVt * conj(Yt(il,:) * diagVnorm) + conj(diagIt) * sparse(1:nl1, t1, Vnorm(t1), nl1, nb);

Sf = V(f1) .* conj(If);
St = V(t1) .* conj(It);
nl_Sf = length(Sf);
 
dAf_dPf = sparse(1:nl_Sf, 1:nl_Sf, 2 * real(Sf), nl_Sf, nl_Sf);
dAf_dQf = sparse(1:nl_Sf, 1:nl_Sf, 2 * imag(Sf), nl_Sf, nl_Sf);
dAt_dPt = sparse(1:nl_Sf, 1:nl_Sf, 2 * real(St), nl_Sf, nl_Sf);
dAt_dQt = sparse(1:nl_Sf, 1:nl_Sf, 2 * imag(St), nl_Sf, nl_Sf);

dAf_dVm = dAf_dPf * real(dSf_dVm) + dAf_dQf * imag(dSf_dVm);
dAf_dVa = dAf_dPf * real(dSf_dVa) + dAf_dQf * imag(dSf_dVa);
dAt_dVm = dAt_dPt * real(dSt_dVm) + dAt_dQt * imag(dSt_dVm);
dAt_dVa = dAt_dPt * real(dSt_dVa) + dAt_dQt * imag(dSt_dVa);

dhn = sparse(2*nl2, n_var);
dhn(:, [1:5 6:10]) = [
    dAf_dVa, dAf_dVm;                     %% "from" flow limit
    dAt_dVa, dAt_dVm;                     %% "to" flow limit
];
%% 
AA = speye(length(x));
ieq = find( abs(xmax-xmin) <= eps );        %% equality constraints
igt = find( xmax >=  1e10 & xmin > -1e10 );     %% greater than, unbounded above
ilt = find( xmin <= -1e10 & xmax <  1e10 );     %% less than, unbounded below
ibx = find( (abs(xmax-xmin) > eps) & (xmax < 1e10) & (xmin > -1e10) );
Ae = AA(ieq, :);
be = xmax(ieq, 1);
Ai  = [ AA(ilt, :); -AA(igt, :); AA(ibx, :); -AA(ibx, :) ];
bi  = [ xmax(ilt, 1); -xmin(igt, 1); xmax(ibx, 1); -xmin(ibx, 1) ];
h = [hn; Ai * x - bi];          %% inequality constraints
g = [gn; Ae * x - be];          %% equality constraints
dh = [dhn' Ai'];                 %% 1st derivative of inequalities
dg = [dgn Ae'];                 %% 1st derivative of equalities
end

%%

 function [f, df] =f_fcn1(x,gen_cost,baseMVA)

[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;
f = sum(x(11:15).*gen_cost(:,COST)*baseMVA); % Finding value of objective function
df = [zeros(10,1);gen_cost(:,COST)/baseMVA;zeros(5,1)]; % first derivatives of objective function
 end
%%


function y = hess_fcn1(x, lam,mu,nb, Ybus,Yf,Yt,Cf,Ct,Sf,St,d2f,dSf_dVa,dSf_dVm,dSt_dVm,dSt_dVa)

Vm = x(6:10,:);
Va = x(1:5,:);

V = Vm .* exp(1j * Va);
%Finding Gxx

lamP = lam(1:5);
lamQ = lam(6:10);
nxtra = length(x) - 2*nb;
[Gpaa, Gpav, Gpva, Gpvv] = d2Sbus_dV2(Ybus, V, lamP);
[Gqaa, Gqav, Gqva, Gqvv] = d2Sbus_dV2(Ybus, V, lamQ);
d2G = [
    real([Gpaa Gpav; Gpva Gpvv]) + imag([Gqaa Gqav; Gqva Gqvv]) sparse(2*nb, nxtra);
    sparse(nxtra, 2*nb + nxtra)
];

% Finding Hxx
muF = mu(1:2);
muT = mu(3:4);

il = [1,6];


[Hfaa, Hfav, Hfva, Hfvv] = d2ASbr_dV2(dSf_dVa, dSf_dVm, Sf, Cf, Yf(il,:), V, muF);
[Htaa, Htav, Htva, Htvv] = d2ASbr_dV2(dSt_dVa, dSt_dVm, St, Ct, Yt(il,:), V, muT);

d2H = [
    [Hfaa Hfav; Hfva Hfvv] + [Htaa Htav; Htva Htvv] sparse(2*nb, nxtra);
    sparse(nxtra, 2*nb + nxtra)
];

y = d2f + d2G + d2H;

end




