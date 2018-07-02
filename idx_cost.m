function [PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost

%% define cost models
PW_LINEAR   = 1;
POLYNOMIAL  = 2;

%% define the indices
MODEL       = 1;    %% cost model, 1 = piecewise linear, 2 = polynomial 
STARTUP     = 2;    %% startup cost in US dollars
SHUTDOWN    = 3;    %% shutdown cost in US dollars
NCOST       = 4;    %% number breakpoints in piecewise linear cost function,
                    %% or number of coefficients in polynomial cost function
COST        = 5;    %% parameters defining total cost function begin in this col
                    
