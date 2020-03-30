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

%% Constraint 55
U55 = zeros(len.Xg,len.U);
for i = 1:len.Xg
    U55(i,data.gen(i,2)) = 1;
end
Aeq55 = zeros(len.Xg,len.total);
Aeq55(:,inp.U) = U55;
Aeq55 = Aeq55(data.statgen == 1,:);
beq55 = ones(size(Aeq55,1),1) .* (1.05^2);

equ(55).Aeq = concA(steps,Aeq55);
equ(55).beq = concB(steps,beq55);

%% Constraint 56
Xbr56 = eye(len.Xbr);
Aeq56 = zeros(len.Xbr,len.total);
Aeq56(:,inp.Xbr) = Xbr56;
Aeq56 = Aeq56(data.statbr == 2,:);
beq56 = zeros(size(Aeq56,1),1);

equ(56).Aeq = initials(steps,Aeq56);
equ(56).beq = beq56;

%% Constraint 57
Xbr57 = eye(len.Xbr);
Aeq57 = zeros(len.Xbr,len.total);
Aeq57(:,inp.Xbr) = Xbr57;
Aeq57 = Aeq57(data.statbr == 0,:);
beq57 = zeros(size(Aeq57,1),1);

equ(57).Aeq = concA(steps,Aeq57);
equ(57).beq = concB(steps,beq57);

%% Constraint 61
Xl61 = eye(len.Xl);
Aeq61 = zeros(len.Xl,len.total);
Aeq61(:,inp.Xl) = Xl61;
Aeq61 = Aeq61(data.statload == 0,:);
beq61 = zeros(size(Aeq61,1),1);

equ(61).Aeq = concA(steps,Aeq61);
equ(61).beq = concB(steps,beq61);

%% Back to normal constraints *SmileyFace
%% Constraint 1
% Searching for bus connection
buscon1 = cell(data.num_bus,6);
for i = 1:data.num_bus
    buscon1{i,1} = data.bus(i,1);
    buscon1{i,2} = data.load(data.load(:,2)==i,1);
    buscon1{i,3} = data.branch(data.branch(:,3)==i,1);
    buscon1{i,4} = data.branch(data.branch(:,2)==i,1);
    buscon1{i,5} = data.gen(data.gen(:,2)==i,1);
    buscon1{i,6} = data.ess(data.ess(:,2)==i,1);
end

Pl1 = zeros(data.num_bus,data.num_load);
Pg1 = zeros(data.num_bus,data.num_gen);
Pbr1 = zeros(data.num_bus,data.num_branch);
Pessc1 = zeros(data.num_bus,data.num_ess);
Pessd1 = Pessc1;

% Input variable constant (-1 and 1) according to the signs of constraint
for i = 1:data.num_bus
    Pl1(buscon1{i,1},buscon1{i,2}) = 1; % P_load,t
    Pbr1(buscon1{i,1},buscon1{i,3}) = -1; % P_hi,t
    Pbr1(buscon1{i,1},buscon1{i,4}) = 1; % P_ij,t
    Pg1(buscon1{i,1},buscon1{i,5}) = -1; % P_g,t
    Pessc1(buscon1{i,1},buscon1{i,6}) = 1; % P_essc,t
    Pessd1(buscon1{i,1},buscon1{i,6}) = -1; % P_essd,t
end
Aeq1 = zeros(len.Sn,len.total);
Aeq1(:,[inp.Pl inp.Pbr inp.Pg inp.Pessc inp.Pessd]) = [Pl1 Pbr1 Pg1 Pessc1 Pessd1];
beq1 = zeros(len.Sn,1);
equ(1).Aeq = concA(steps,Aeq1);
equ(1).beq = concB(steps,beq1);

%% Constraint 2
Aeq2 = zeros(len.Sn,len.total);
Aeq2(:,[inp.Ql inp.Qbr inp.Qg inp.Qessc inp.Qessd]) = [Pl1 Pbr1 Pg1 Pessc1 Pessd1]; % Same construction
beq2 = zeros(len.Sn,1);
equ(2).Aeq = concA(steps,Aeq2);
equ(2).beq = concB(steps,beq2);

%% Running the MILP
RunMILP;
toc;