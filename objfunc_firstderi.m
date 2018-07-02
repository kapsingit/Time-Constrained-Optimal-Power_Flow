
function [f,df] = objfunc_firstderi(x, T)

[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

[bus_data,gen_data,branch_data, gen_cost,baseMVA] = network_info;

f_value = 0;

for i = 1:T
    f_value = f_value +sum(x(11*i:11*i+4).*gen_cost(:,COST)*baseMVA); % Finding value of objective function
end
f = f_value; 
df = [zeros(10,1);gen_cost(:,COST)/baseMVA;zeros(5,1)]; % first derivatives of objective function
df = repmat(df,T,1);
