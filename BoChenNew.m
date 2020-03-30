%% MILP Service Restoration BoChen
tic; 
clear;
steps = 10; % maximum time step
intmin = 1; % in minute
inthour = intmin/60; % in hour
bigM = 20000;

%% Retrieve data
data = data13bochenrev;

%% Distflow equation
len = lengthbochen(data);
inp = inputvariables(len);

%% INITIAL CONDITION CONSTRAINTS
%% Constraint 54 // can be used as bound
Xg54 = eye(len.Xg);
Aeq54 = zeros(len.Xg,len.total);
Aeq54(:,inp.Xg) = Xg54;
Aeq54 = Aeq54(data.statgen == 1,:);
beq54 = ones(size(Aeq54,1),1);

equ(54).Aeq = concA(steps,Aeq54);
equ(54).beq = concB(steps,beq54);
