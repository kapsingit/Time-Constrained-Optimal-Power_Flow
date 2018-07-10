
%% Code generation
%% kapil.duwadi@jacks.sdstate.edu
%% Time constrained optimal power flow (TCOPF)

% A 5 bus PJM transmission network is being used as a benchmark 
% All data gathered from matpower ("case5.m")

%% 7:52 AM 7/1/2018 # Starting attempt
clear;clc;
% Reads the bus, branch, gen and cost information about the network

%% Gathering data for PJM 5 bus system

[bus_data,gen_data,branch_data, gen_cost,baseMVA] = network_info;

%%  Making the addmittance network 

[Ybus, Yf, Yt] = make_ybus;

%% Initiallizing the decision variables

timesteps = 24;
[x, xmin, xmax] = decesion_variables(timesteps);

%% Finding objective function and it's first derivative

%First stop: 8:51 AM 7/1/2018
% Resumes at 12:18 PM 7/1/2018

[f,df] = objfunc_firstderi(x,timesteps);

%% Finding second derivative of objective function

d2f = sparse(1:20*timesteps,1:20*timesteps,0); % Using linear cost functions

%% Finding equality and ineqiality constraints and their first derivatives

% Equality constraints 
        % 1) Power mismatch at node is equal to zero
        % 2) voltage angle at reference node is equal to zero

        % Number of equality constraints = 2T*nb+T = 2*24*5+24 = 264 (size
        % of g)
% Ineqiality constraimtts

        % 1) line limits
        % 2) decision variables must be greater than lower bound and lesser
        % than upper bound
        
        % Number of inequality constraints = 2*T*nl+6*nb*T = 2*24*2+6*5*24
        % = 816 (size of h) Note: No limit on voltage angles also only two
        % lines have limit in this case
% second stop: 12:34 PM 7/1/2018
% resumes at 3:57 PM 7/1/2018

[load_data_p, load_data_q] = generating_load_data(timesteps);

[g,gn,h,hn,dh,dg,Sflist,Stlist,dSf_dValist,dSf_dVmlist,dSt_dVmlist,dSt_dValist] = Eqcons_df(x,xmax,xmin,timesteps,load_data_p,load_data_q);

%% Grabing some dimensions and initiallizing constants of the optimization
neq = size(g, 1);           %% number of equality constraints
niq = size(h, 1);           %% number of inequality constraints
neqnln = size(gn, 1);       %% number of nonlinear equality constraints
niqnln = size(hn, 1);       %% number of nonlinear inequality constraints

% initialize gamma, lam, mu, z, e
z0 = 1;
gamma = 1;                  
lam = zeros(neq, 1);
z   = z0 * ones(niq, 1);
mu  = z;
k = find(h < -z0);
z(k) = -h(k);
e = ones(niq, 1);
%% Solving optimization
kk = 1;
err = 1;
cond_TCOPF = [];
while kk <= 500
    zinvdiag = sparse(1:niq, 1:niq, 1./ z, niq, niq);
    mudiag = sparse(1:niq, 1:niq, mu, niq, niq);
    dh_zinv = dh * zinvdiag;
    Lx = df + dg * lam + dh * mu;
    Lxx = hess_func(x,lam,mu,timesteps,d2f,Sflist,Stlist,dSf_dValist,dSf_dVmlist,dSt_dVmlist,dSt_dValist);
    M = Lxx + dh_zinv * mudiag * dh';
    N = Lx + dh_zinv * (mudiag * h + gamma * e);
    W = [M dg;dg' sparse(neq, neq)];
    B = [-N; -g];

    dxdlam = W \ B;
    
    nx = size(x,1);
    dx = dxdlam(1:nx, 1);
    dlam = dxdlam(nx+(1:neq), 1);
    dz = -h - z - dh' * dx;
    dmu = -mu + zinvdiag *(gamma*e - mudiag * dz);
    
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
        alphad = min( [xi * min(mu(k) ./ -dmu(k)) 1] );
    end
    
    x = x + alphap * dx;
    z = z + alphap * dz;
    lam = lam + alphad * dlam;
    mu  = mu  + alphad * dmu;
    
    
    if niq > 0
        gamma = 0.1* (z' * mu) / niq;
    end
    [f,df] = objfunc_firstderi(x,timesteps);              %% cost
    [g,gn,h,hn,dh,dg,Sflist,Stlist,dSf_dValist,dSf_dVmlist,dSt_dVmlist,dSt_dValist] = Eqcons_df(x,xmax,xmin,timesteps,load_data_p,load_data_q);

   solu(kk,:) = x;
   if kk>2
   err = max(solu(kk,:)-solu(kk-1,:));
   end
   if err<0.000001
       break;
   end
   kk = kk+1;
   cond_TCOPF = [cond_TCOPF cond(full(W))];
end
%% Post processing
sol = reshape(x,20,timesteps);
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
%%
figure('color',[1,1,1]);
plot(solu(2:end,:)-solu(1:end-1,:));
xlabel('Iteration');
ylabel('Error magnitude');
grid on; box on;
%% Error plotting