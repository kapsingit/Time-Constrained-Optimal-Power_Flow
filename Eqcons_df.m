
function [gt,gnt,ht,hnt,dht,dgt,Sflist,Stlist,dSf_dValist,dSf_dVmlist,dSt_dVmlist,dSt_dValist] = Eqcons_df(xt,xmaxt,xmint,T,load_data_p,load_data_q)

[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, RAMP] = idfor_gen;

[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, VM, ...
    VA, BASE_KV, VMAX, VMIN] = idfor_bus;

[bus_data,gen_data,branch_data, gen_cost,baseMVA] = network_info;

[Ybus, Yf, Yt] = make_ybus;

[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C,ANGMIN, ANGMAX ] = idfor_branch;

ng = size(gen_data,1);
nb = size(bus_data,1);

gen_buses = gen_data(:,GEN_BUS);

Cg = sparse(gen_buses, (1:ng)', 1, nb, ng); % generator connection matrix

gnt = [];
hnt = [];
dgnt = [];
dhnt = [];
dht = [];
dgt = [];
gt = [];
ht = [];

Sflist = struct;
Stlist = struct;
dSf_dValist = struct;
dSf_dVmlist = struct;
dSt_dVmlist = struct;
dSt_dValist = struct;

for i = 1:T
    
    x = xt([1:20]+(i-1)*20);
    xmax = xmaxt([1:20]+(i-1)*20);
    xmin = xmint([1:20]+(i-1)*20);
    Va = x(1:5,:);
    Vm = x(6:10,:);
    Pg = x(11:15,:);
    Qg = x(16:20,:);
    
    Sbusg = Cg*(Pg + 1j*Qg);
    Sload = (load_data_p(:,i)+1j*load_data_q(:,i))/baseMVA;
    Sbus = Sbusg - Sload;
    
    V = Vm .* exp(1j * Va);
    
    mis = V .* conj(Ybus * V) - Sbus;
    
    gn = [ real(mis);            %% active power mismatch for all buses
      imag(mis) ];          %% reactive power mismatch for all buses

    gnt = [gnt;gn];
    
    il = find(branch_data(:,RATE_A)~=0);
    flow_max = (branch_data(il, RATE_A)/baseMVA).^2;
    nl2 = length(il);
    
    Sf = V(branch_data(il, F_BUS)) .* conj(Yf(il,:) * V);
    St = V(branch_data(il, T_BUS)) .* conj(Yt(il,:) * V);
    
    hn = [ Sf .* conj(Sf) - flow_max;      %% branch apparent power limits (from bus)
      St .* conj(St) - flow_max ];
    hnt = [hnt; hn];
  
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
  dgnt = blkdiag(dgnt,dgn); % changed
  dhnt = blkdiag(dhnt,dhn); % changed
  
  Sflist(i).value = Sf;
  Stlist(i).value = St;
  dSf_dValist(i).value = dSf_dVa;
  dSf_dVmlist(i).value = dSf_dVm;
  dSt_dVmlist(i).value = dSt_dVm;
  dSt_dValist(i).value = dSt_dVa;
  
  AA = speye(length(x));
  ieq = find( abs(xmax-xmin) <= eps );        %% equality constraints
  igt = find( xmax >=  1e10 & xmin > -1e10 );     %% greater than, unbounded above
  ilt = find( xmin <= -1e10 & xmax <  1e10 );     %% less than, unbounded below
  ibx = find( (abs(xmax-xmin) > eps) & (xmax < 1e10) & (xmin > -1e10) );
  Ae = AA(ieq, :);
  be = xmax(ieq, 1);
  Ai  = [ AA(ilt, :); -AA(igt, :); AA(ibx, :); -AA(ibx, :) ];
  bi  = [ xmax(ilt, 1); -xmin(igt, 1); xmax(ibx, 1); -xmin(ibx, 1) ];
  ht_a = [hn; Ai * x - bi];          %% inequality constraints
  gt_a = [gn; Ae * x - be];          %% equality constraints
  ht = [ht;ht_a];
  gt = [gt;gt_a];
  dht_a = [dhn' Ai'];                 %% 1st derivative of inequalities
  dgt_a = [dgn Ae'];                 %% 1st derivative of equalities
  dht = blkdiag(dht,dht_a); 
  dgt = blkdiag(dgt,dgt_a);
end

% Limits

%   AA = speye(length(xt));
%   ieq = find( abs(xmaxt-xmint) <= eps );        %% equality constraints
%   igt = find( xmaxt >=  1e10 & xmint > -1e10 );     %% greater than, unbounded above
%   ilt = find( xmint <= -1e10 & xmaxt <  1e10 );     %% less than, unbounded below
%   ibx = find( (abs(xmaxt-xmint) > eps) & (xmaxt < 1e10) & (xmint > -1e10) );
%   Ae = AA(ieq, :);
%   be = xmaxt(ieq, 1);
%   Ai  = [ AA(ilt, :); -AA(igt, :); AA(ibx, :); -AA(ibx, :) ];
%   bi  = [ xmaxt(ilt, 1); -xmint(igt, 1); xmaxt(ibx, 1); -xmint(ibx, 1) ];
%   ht = [hnt;Ai * xt - bi];
%   gt = [gnt; Ae * xt - be];
%   dht = [dhnt' Ai']; 
%   dgt = [dgnt Ae'];


% Ramping constraints and its derivative
id = [];
for i = 11:15
    for j = 1:T-1
        id = [id;i+(j-1)*20 i+j*20]; % 20 is for 5 bus is actually  (total number of decision variables)
    end
end
num_rows = (T-1)*ng;
Afirst = sparse(repmat([1:num_rows],1,2),id,[ones(num_rows,1) -ones(num_rows,1)],num_rows,length(xt));
Asecond = sparse(repmat([1:num_rows],1,2),id,[-ones(num_rows,1) ones(num_rows,1)],num_rows,length(xt));
bfs = repmat(gen_data(:,RAMP)/baseMVA,1,T-1);
bfs = bfs(:);
ht = [ht;Afirst*xt-bfs;Asecond*xt-bfs];
dht = [dht Afirst' Asecond'];
