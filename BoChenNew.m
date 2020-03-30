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

%% Constraint 52
Xbr51 = eye(len.Xbr);
Sn51_bef = zeros(len.Xbr,len.Sn);
for i = 1:len.Xbr
    Sn51_bef(i,data.branch(i,2)) = -1;
    Sn51_bef(i,data.branch(i,3)) = -1;
end
A52 = zeros(len.Xbr,len.total);
A52_bef = A52;
A52(:,inp.Xbr) = Xbr51;
A52 = A52(data.statbr == 2,:);

A52_bef(:,inp.Sn) = Sn51_bef;
A52_bef = A52_bef(data.statbr == 2,:);

b52 = zeros(size(A52,1),1);
[ineq(52).A,ineq(52).b] = time_relate(steps,A52_bef,A52,b52);

%% Constraint 53
Sn53 = ones(1,len.Sn);
Xbr53 = -ones(1,len.Xbr);
Xg53 = -ones(1,len.Xg);
Xg53(data.statgen ~= 1) = 0;

Aeq53 = zeros(1,len.total);
Aeq53(:,inp.Sn) = Sn53;
Aeq53(:,inp.Xbr) = Xbr53;
Aeq53(:,inp.Xg) = Xg53;

beq53 = zeros(size(Aeq53,1));

equ(53).Aeq = concA(steps,Aeq53);
equ(53).beq = concB(steps,beq53);

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

%% Constraint 3
% C.3 -- A (( -M.X_br - P_br <= 0 ))
Xbr3a = -bigM * eye(len.Xbr);
Pbr3a = -eye(len.Pbr);
A3a = zeros(size(Pbr3a,1),len.total);
A3a(:,inp.Pbr) = Pbr3a; 
A3a(:,inp.Xbr) = Xbr3a;
b3a = zeros(size(Pbr3a,1),1);

% C.3 -- B (( P_br - M.X_br <= 0 ))
Xbr3b = -bigM * eye(len.Xbr);
Pbr3b = -Pbr3a;
A3b = zeros(size(Pbr3b,1),len.total);
A3b(:,inp.Pbr) = Pbr3b; 
A3b(:,inp.Xbr) = Xbr3b;
b3b = b3a;

% Concatenate
A3 = [A3a; A3b]; b3 = [b3a; b3b];
ineq(3).A = concA(steps,A3);
ineq(3).b = concB(steps,b3);

%% Constraint 4
% C.3 -- A (( -M.X_br - Q_br <= 0 ))
Xbr4a = -bigM * eye(len.Xbr);
Qbr4a = -eye(len.Qbr);
A4a = zeros(size(Qbr4a,1),len.total);
A4a(:,inp.Qbr) = Qbr4a;
A4a(:,inp.Xbr) = Xbr4a;
b4a = zeros(size(Qbr4a,1),1);

% C.3 -- B (( Q_br - M.X_br <= 0 ))
Xbr4b = -bigM * eye(len.Xbr);
Qbr4b = eye(len.Qbr);
A4b = zeros(size(Qbr4b,1),len.total);
A4b(:,inp.Qbr) = Qbr4b;
A4b(:,inp.Xbr) = Xbr4b;
b4b = b4a;

A4 = [A4a; A4b];
b4 = [b4a; b4b];
ineq(4).A = concA(steps,A4);
ineq(4).b = concB(steps,b4);


%% Constraint 5 (( U_i - U_j - 2.r_ij.P_br - 2.x_ij.Qbr + M.Xbr <= M ))
U5 = zeros(len.Pbr,len.Sn);
Pbr5 = zeros(len.Pbr); Qbr5 = Pbr5;
Xbr5 = zeros(len.Xbr);
for i = 1:len.Pbr
    U5(i,data.branch(i,2)) = 1; % U_it
    U5(i,data.branch(i,3)) = -1; % U_jt
    Pbr5(i,i) = -2 * data.branch(i,6); % - 2.r_ij.P_brt
    Qbr5(i,i) = -2 * data.branch(i,7); % - 2.x_ij.Qbrt
    Xbr5(i,i) = bigM; % M.Xbrt
end
A5 = zeros(size(U5,1),len.total);
A5(:,(inp.U)) = U5;
A5(:,(inp.Pbr)) = Pbr5;
A5(:,(inp.Qbr)) = Qbr5;
A5(:,(inp.Xbr)) = Xbr5;
b5 = ones(size(U5,1),1) .* bigM;
ineq(5).A = concA(steps,A5);
ineq(5).b = concB(steps,b5);

%% Constraint 6 (( - U_i + U_j + 2.r_ij.P_br + 2.x_ij.Qbr + M.Xbr <= M ))
% Be careful, Xbr does not change
U6 = -U5; Pbr6 = -Pbr5; Qbr6 = -Qbr5; Xbr6 = Xbr5;
A6 = zeros(size(U6,1),len.total);
A6(:,inp.U) = U6;
A6(:,inp.Pbr) = Pbr6;
A6(:,inp.Qbr) = Qbr6;
A6(:,inp.Xbr) = Xbr6;
b6 = b5;
ineq(6).A = concA(steps,A6);
ineq(6).b = concB(steps,b6);

%% Constraint 7-8 only for CLPU, skipped
%% Constraint 9 and 10
Pl9 = eye(len.Pl); Xl9 = -eye(len.Xl) .* data.load(:,4);
Ql10 = eye(len.Ql); Xl10 = -eye(len.Xl) .* data.load(:,5);

Aeq9 = zeros(size(Pl9,1),len.total);
Aeq9(:,inp.Pl) = Pl9;
Aeq9(:,inp.Xl) = Xl9;

Aeq10 = zeros(size(Ql10,1),len.total);
Aeq10(:,inp.Ql) = Ql10;
Aeq10(:,inp.Xl) = Xl10;

beq9 = zeros(len.Xl,1);
beq10 = zeros(len.Xl,1);

equ(9).Aeq = concA(steps,Aeq9);
equ(9).beq = concB(steps,beq9);
equ(10).Aeq = concA(steps,Aeq10);
equ(10).beq = concB(steps,beq10);

%% Constraint 13
Sijmax = data.branch(:,9);
Sij = Sijmax.*sqrt((2*pi/6)/sin(2*pi/6));

% C.13 -- A 
Pbr13a = -sqrt(3) * eye(len.Pbr);
Qbr13a = -eye(len.Qbr);
A13a = zeros(size(Pbr13a,1),len.total);
A13a(:,inp.Pbr) = Pbr13a;
A13a(:,inp.Qbr) = Qbr13a;
b13a = sqrt(3) .* Sij;

% C.13 -- B
Pbr13b = -Pbr13a; Qbr13b = -Qbr13a;
A13b = zeros(size(Pbr13b,1),len.total);
A13b(:,inp.Pbr) = Pbr13b;
A13b(:,inp.Qbr) = Qbr13b;
b13b = b13a;

A13 = [A13a; A13b];
b13 = [b13a; b13b];

ineq(13).A = concA(steps,A13);
ineq(13).b = concB(steps,b13);

%% Constraint 14
% C.14 -- A
Qbr14a = -eye(len.Qbr);
A14a = zeros(size(Qbr14a,1),len.total);
A14a(:,(inp.Qbr)) = Qbr14a;
b14a = b13a/2;

% C.14 -- B
Qbr14b = eye(len.Qbr);
A14b = zeros(size(Qbr14b,1),len.total);
A14b(:,(inp.Qbr)) = Qbr14b;
b14b = b14a;

A14 = [A14a; A14b];
b14 = [b14a; b14b];

ineq(14).A = concA(steps,A14);
ineq(14).b = concB(steps,b14);


%% Running the MILP
RunMILP;
toc;