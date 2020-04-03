%% MILP Service Restoration BoChen
tic; 
close all;
clear;
setVmax = 1.05; % set maximum voltage
setVmin = 0.95; % set minimum voltage
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
Xg53(data.statgen ~= 1) = 0; % Set the non black-start DG to 0

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
beq55 = ones(size(Aeq55,1),1) .* (setVmax^2);

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
A14a(:,inp.Qbr) = Qbr14a;
b14a = b13a/2;

% C.14 -- B
Qbr14b = eye(len.Qbr);
A14b = zeros(size(Qbr14b,1),len.total);
A14b(:,inp.Qbr) = Qbr14b;
b14b = b14a;

A14 = [A14a; A14b];
b14 = [b14a; b14b];

ineq(14).A = concA(steps,A14);
ineq(14).b = concB(steps,b14);

%% Constraint 15
% C.15 -- A
Pbr15a = sqrt(3) * eye(len.Pbr);
Qbr15a = -eye(len.Qbr);
A15a = zeros(size(Pbr15a,1),len.total);
A15a(:,inp.Pbr) = Pbr15a;
A15a(:,inp.Qbr) = Qbr15a;
b15a = sqrt(3) .* Sij;

% C.15 -- B
Pbr15b = -Pbr15a; Qbr15b = -Qbr15a;
A15b = zeros(size(Pbr15b,1),len.total);
A15b(:,inp.Pbr) = Pbr15b;
A15b(:,inp.Qbr) = Qbr15b;
b15b = b15a;

A15 = [A15a; A15b];
b15 = [b15a; b15b];

ineq(15).A = concA(steps,A15);
ineq(15).b = concB(steps,b15);

%% Constraint 16
spin = 15/100; % Spinning reserve is set to 15 percent
Pl16 = (1+spin) .* ones(1,len.Pl);
Xg16 = -data.gen(:,6)';
Xessd16 = -data.ess(:,15)';
A16 = zeros(size(Pl16,1),len.total);
A16(:,inp.Pl) = Pl16;
A16(:,inp.Xg) = Xg16;
A16(:,inp.Xessd) = Xessd16;
b16 = 0;

ineq(16).A = concA(steps,A16);
ineq(16).b = concB(steps,b16);

%% Constraint 17
% C.17 -- A
Pg17a = -eye(len.Pg);
Xg17a = eye(len.Xg) .* data.gen(:,7);
A17a = zeros(size(Pg17a,1),len.total);
A17a(:,inp.Pg) = Pg17a;
A17a(:,inp.Xg) = Xg17a;
b17a = zeros(len.Xg,1);

% C.17 -- B
Pg17b = eye(len.Pg);
Xg17b = -eye(len.Xg) .* data.gen(:,6);
A17b = zeros(size(Pg17b,1),len.total);
A17b(:,inp.Pg) = Pg17b;
A17b(:,inp.Xg) = Xg17b;
b17b = b17a;

A17 = [A17a; A17b];
b17 = [b17a; b17b];

ineq(17).A = concA(steps,A17);
ineq(17).b = concB(steps,b17);

%% Constraint 18
% C.18 -- A
Qg18a = -eye(len.Pg);
Xg18a = eye(len.Xg) .* data.gen(:,9);
A18a = zeros(size(Qg18a,1),len.total);
A18a(:,inp.Qg) = Qg18a;
A18a(:,inp.Xg) = Xg18a;
b18a = b17a;

% C.18 -- B
Qg18b = eye(len.Pg);
Xg18b = -eye(len.Xg) .* data.gen(:,8);
A18b = zeros(size(Qg18b,1),len.total);
A18b(:,inp.Qg) = Qg18b;
A18b(:,inp.Xg) = Xg18b;
b18b = b18a;

A18 = [A18a; A18b];
b18 = [b18a; b18b];

ineq(18).A = concA(steps,A18);
ineq(18).b = concB(steps,b18);

%% Constraint 19
Pg19 = eye(len.Pg) .* tan(acos(data.gen(:,4)));
Qg19 = -eye(len.Qg);

Aeq19 = zeros(size(Pg19,1),len.total);
Aeq19(:,inp.Pg) = Pg19;
Aeq19(:,inp.Qg) = Qg19;
Aeq19 = Aeq19(data.statgen ~= 1,:);
beq19 = zeros(size(Aeq19,1),1);

equ(19).Aeq = concA(steps,Aeq19);
equ(19).beq = concB(steps,beq19);

%% Constraint 20
% C.20 -- A
Pessc20a = -eye(len.Pessc);
Xessc20a = eye(len.Xessc) .* data.ess(:,10);
A20a = zeros(len.Pessc,len.total);
A20a(:,inp.Pessc) = Pessc20a;
A20a(:,inp.Xessc) = Xessc20a;
b20a = zeros(len.Pessc,1);

% C.20 -- B
Pessc20b = eye(len.Pessc);
Xessc20b = -eye(len.Xessc) .* data.ess(:,11);
A20b = zeros(len.Pessc,len.total);
A20b(:,inp.Pessc) = Pessc20b;
A20b(:,inp.Xessc) = Xessc20b;
b20b = b20a;

A20 = [A20a; A20b];
b20 = [b20a; b20b];

ineq(20).A = concA(steps,A20);
ineq(20).b = concB(steps,b20);

%% Constraint 21
% C.21 -- A
Pessd21a = -eye(len.Pessd);
Xessd21a = eye(len.Xessd) .* data.ess(:,14);
A21a = zeros(len.Pessd,len.total);
A21a(:,inp.Pessd) = Pessd21a;
A21a(:,inp.Xessd) = Xessd21a;
b21a = zeros(len.Pessd,1);

% C.21 -- B
Pessd21b = eye(len.Pessd);
Xessd21b = -eye(len.Xessd) .* data.ess(:,15);
A21b = zeros(len.Pessd,len.total);
A21b(:,inp.Pessd) = Pessd21b;
A21b(:,inp.Xessd) = Xessd21b;
b21b = zeros(len.Pessd,1);

A21 = [A21a; A21b];
b21 = [b21a; b21b];

ineq(21).A = concA(steps,A21);
ineq(21).b = concB(steps,b21);

%% Constraint 22
% C.22 -- A
Qessc22a = -eye(len.Qessc);
Xessc22a = eye(len.Xessc) .* data.ess(:,12);
A22a = zeros(len.Qessc,len.total);
A22a(:,(inp.Qessc)) = Qessc22a;
A22a(:,(inp.Xessc)) = Xessc22a;
b22a = zeros(len.Qessc,1);

% C.22 -- B
Qessc22b = eye(len.Qessc);
Xessc22b = -eye(len.Xessc) .* data.ess(:,13);
A22b = zeros(len.Qessc,len.total);
A22b(:,(inp.Qessc)) = Qessc22b;
A22b(:,(inp.Xessc)) = Xessc22b;
b22b = zeros(len.Qessc,1);

A22 = [A22a; A22b];
b22 = [b22a; b22b];

ineq(22).A = concA(steps,A22);
ineq(22).b = concB(steps,b22);

%% Constraint 23
% C.23 -- A
Qessd23a = -eye(len.Qessd);
Xessd23a = eye(len.Xessd) .* data.ess(:,16);
A23a = zeros(len.Qessd,len.total);
A23a(:,(inp.Qessd)) = Qessd23a;
A23a(:,(inp.Xessd)) = Xessd23a;
b23a = zeros(len.Qessd,1);

% C.23 -- B
Qessd23b = eye(len.Qessd);
Xessd23b = -eye(len.Xessd) .* data.ess(:,17);
A23b = zeros(len.Qessd,len.total);
A23b(:,(inp.Qessd)) = Qessd23b;
A23b(:,(inp.Xessd)) = Xessd23b;
b23b = zeros(len.Qessd,1);

A23 = [A23a; A23b];
b23 = [b23a; b23b];

ineq(23).A = concA(steps,A23);
ineq(23).b = concB(steps,b23);

%% Constraint 24
Xessc24 = eye(len.Xessc);
Xessd24 = Xessc24;
Sn24 = zeros(len.Xessc,len.Sn);
for i = 1:len.Xessc
    Sn24(i,data.ess(i,2)) = -1;
end
A24 = zeros(len.Xessc,len.total);
A24(:,inp.Sn) = Sn24;
A24(:,inp.Xessc) = Xessc24;
A24(:,inp.Xessd) = Xessd24;
b24 = zeros(len.Xessc,1);

ineq(24).A = concA(steps,A24);
ineq(24).b = concB(steps,b24);

%% Constraint 25 // Time dependent, initial condition
Eess25 = eye(len.Eess);
Aeq25 = zeros(len.Eess,len.total);
Aeq25(:,inp.Eess) = Eess25;
beq25 = data.ess(:,5).*data.ess(:,4);

equ(25).Aeq = initials(steps,Aeq25);
equ(25).beq = beq25;

%% Constraint 26 // Time dependent
Eess26 = eye(len.Eess);
Pessc26 = -inthour*eye(len.Pessc).*data.ess(:,8);
Pessd26 = inthour*eye(len.Pessd)./data.ess(:,9);
Aeq26 = zeros(len.Eess,len.total);
Aeq26_bef = Aeq26;
Aeq26(:,inp.Eess) = Eess26;
Aeq26(:,inp.Pessc) = Pessc26;
Aeq26(:,inp.Pessd) = Pessd26;

Eess26_bef = -eye(len.Eess);
Aeq26_bef(:,inp.Eess) = Eess26_bef;

beq26 = zeros(len.Eess,1);
[equ(26).Aeq, equ(26).beq] = time_relate(steps,Aeq26_bef,Aeq26,beq26);

%% Constraint 27 // Can be changed directly into bound
% C.27 -- A
Eess27a = -eye(len.Eess);
A27a = zeros(len.Eess,len.total);
A27a(:,inp.Eess) = Eess27a;
b27a = -data.ess(:,4) .* data.ess(:,6);

% C.27 -- B
Eess27b = eye(len.Eess);
A27b = zeros(len.Eess,len.total);
A27b(:,inp.Eess) = Eess27b;
b27b = data.ess(:,4) .* data.ess(:,7);

A27 = [A27a; A27b];
b27 = [b27a; b27b];

ineq(27).A = concA(steps,A27);
ineq(27).b = concB(steps,b27);

%% Constraint 28
Pessc28 = eye(len.Pessc);
A28 = zeros(len.Pessc,len.total);
A28(:,inp.Pessc) = Pessc28;
b28 = data.ess(:,18)*intmin;

ineq(28).A = initials(steps,A28);
ineq(28).b = b28;

%% Constraint 29
Pessd29 = eye(len.Pessd);
A29 = zeros(len.Pessd,len.total);
A29(:,inp.Pessd) = Pessd29;
b29 = data.ess(:,19)*intmin;

ineq(29).A = initials(steps,A29);
ineq(29).b = b29;

%% Constraint 30
Qessc30 = eye(len.Qessc);
A30 = zeros(len.Qessc,len.total);
A30(:,inp.Qessc) = Qessc30;
b30 = data.ess(:,20)*intmin;

ineq(30).A = initials(steps,A30);
ineq(30).b = b30;

%% Constraint 31
Qessd31 = eye(len.Qessd);
A31 = zeros(len.Qessd,len.total);
A31(:,inp.Qessd) = Qessd31;
b31 = data.ess(:,21)*intmin;

ineq(31).A = initials(steps,A31);
ineq(31).b = b31;

%% Constraint 32
% C.32 -- A
Pessc32a = -eye(len.Pessc);
Pessc32a_bef = eye(len.Pessc);
A32a = zeros(len.Pessc,len.total);
A32a_bef = A32a;
A32a(:,inp.Pessc) = Pessc32a;
A32a_bef(:,inp.Pessc) = Pessc32a_bef;
b32a = intmin .* data.ess(:,18);

[A32a, b32a] = time_relate(steps,A32a_bef,A32a,b32a);

% C.32 -- B
Pessc32b = eye(len.Pessc);
Pessc32b_bef = -eye(len.Pessc);
A32b = zeros(len.Pessc,len.total);
A32b_bef = A32b;
A32b(:,inp.Pessc) = Pessc32b;
A32b_bef(:,inp.Pessc) = Pessc32b_bef;
b32b = intmin .* data.ess(:,18);

[A32b, b32b] = time_relate(steps,A32b_bef,A32b,b32b);

A32 = [A32a; A32b];
b32 = [b32a; b32b];

ineq(32).A = A32;
ineq(32).b = b32;

%% Constraint 33
% C.33 -- A
Pessd33a = -eye(len.Pessd);
Pessd33a_bef = eye(len.Pessd);
A33a = zeros(len.Pessd,len.total);
A33a_bef = A33a;
A33a(:,inp.Pessd) = Pessd33a;
A33a_bef(:,inp.Pessd) = Pessd33a_bef;
b33a = intmin .* data.ess(:,19);

[A33a, b33a] = time_relate(steps,A33a_bef,A33a,b33a);

% C.33 -- B
Pessd33b = eye(len.Pessd);
Pessd33b_bef = -eye(len.Pessd);
A33b = zeros(len.Pessd,len.total);
A33b_bef = A33b;
A33b(:,inp.Pessd) = Pessd33b;
A33b_bef(:,inp.Pessd) = Pessd33b_bef;
b33b = intmin .* data.ess(:,19);

[A33b, b33b] = time_relate(steps,A33b_bef,A33b,b33b);

A33 = [A33a; A33b];
b33 = [b33a; b33b];

ineq(33).A = A33;
ineq(33).b = b33;

%% Constraint 34
% C.34 -- A
Qessc34a = -eye(len.Qessc);
Qessc34a_bef = eye(len.Qessc);
A34a = zeros(len.Qessc,len.total);
A34a_bef = A34a;
A34a(:,inp.Qessc) = Qessc34a;
A34a_bef(:,inp.Qessc) = Qessc34a_bef;
b34a = intmin .* data.ess(:,20);

[A34a, b34a] = time_relate(steps,A34a_bef,A34a,b34a);

% C.34 -- B
Qessc34b = eye(len.Qessc);
Qessc34b_bef = -eye(len.Qessc);
A34b = zeros(len.Qessc,len.total);
A34b_bef = A34b;
A34b(:,inp.Qessc) = Qessc34b;
A34b_bef(:,inp.Qessc) = Qessc34b_bef;
b34b = intmin .* data.ess(:,20);

[A34b, b34b] = time_relate(steps,A34b_bef,A34b,b34b);

A34 = [A34a; A34b];
b34 = [b34a; b34b];

ineq(34).A = A34;
ineq(34).b = b34;

%% Constraint 35
% C.35 -- A
Qessd35a = -eye(len.Qessd);
Qessd35a_bef = eye(len.Qessd);
A35a = zeros(len.Qessd,len.total);
A35a_bef = A35a;
A35a(:,inp.Qessd) = Qessd35a;
A35a_bef(:,inp.Qessd) = Qessd35a_bef;
b35a = intmin .* data.ess(:,21);

[A35a, b35a] = time_relate(steps,A35a_bef,A35a,b35a);

% C.35 -- B
Qessd35b = eye(len.Qessd);
Qessd35b_bef = -eye(len.Qessd);
A35b = zeros(len.Qessd,len.total);
A35b_bef = A35b;
A35b(:,inp.Qessd) = Qessd35b;
A35b_bef(:,inp.Qessd) = Qessd35b_bef;
b35b = intmin .* data.ess(:,21);

[A35b, b35b] = time_relate(steps,A35b_bef,A35b,b35b);

A35 = [A35a; A35b];
b35 = [b35a; b35b];

ineq(35).A = A35;
ineq(35).b = b35;

%% Constraint 36
con36a = zeros(len.Pessc,len.total);
con36b = con36a; con36c = con36a; con36d = con36a;

Pessc36a = eye(len.Pessc); 
Pessd36b = eye(len.Pessd);  
Qessc36c = eye(len.Qessc);  
Qessd36d = eye(len.Qessd); 

con36a(:,inp.Pessc) = Pessc36a;
con36b(:,inp.Pessd) = Pessd36b;
con36c(:,inp.Qessc) = Qessc36c;
con36d(:,inp.Qessd) = Qessd36d;

Aeq36a = lasts(steps,con36a); beq36a = zeros(len.Pessc,1);
Aeq36b = lasts(steps,con36b); beq36b = zeros(len.Pessd,1);
Aeq36c = lasts(steps,con36c); beq36c = zeros(len.Qessc,1);
Aeq36d = lasts(steps,con36d); beq36d = zeros(len.Qessd,1);

Aeq36 = [Aeq36a; Aeq36b; Aeq36c; Aeq36d];
beq36 = [beq36a; beq36b; beq36c; beq36d];

equ(36).Aeq = Aeq36;
equ(36).beq = beq36;

%% Constraint 37 
% C.37 -- A
U37a = -eye(len.U);
Sn37a = eye(len.Sn) .* (setVmin^2);
A37a = zeros(len.U,len.total);
A37a(:,inp.U) = U37a;
A37a(:,inp.Sn) = Sn37a;
b37a = zeros(len.U,1);

% C.37 -- B
U37b = eye(len.U);
Sn37b = -eye(len.Sn) .* (setVmax^2);
A37b = zeros(len.U,len.total);
A37b(:,inp.U) = U37b;
A37b(:,inp.Sn) = Sn37b;
b37b = zeros(len.U,1);

A37 = [A37a; A37b];
b37 = [b37a; b37b];

ineq(37).A = concA(steps,A37);
ineq(37).b = concB(steps,b37);

%% Constraint 38 // time dependent
% C.38 -- A
Pg38a = -eye(len.Pg);
Pg38a_bef = eye(len.Pg);
A38a = zeros(len.Pg,len.total);
A38a(:,inp.Pg) = Pg38a;
A38a_bef = zeros(len.Pg,len.total);
A38a_bef(:,inp.Pg) = Pg38a_bef;
b38a = data.gen(:,10) * intmin;
[A38a,b38a] = time_relate(steps,A38a_bef,A38a,b38a);

% C.38 -- B
Pg38b = eye(len.Pg);
Pg38b_bef = -eye(len.Pg);
A38b = zeros(len.Pg,len.total);
A38b(:,inp.Pg) = Pg38b;
A38b_bef = zeros(len.Pg,len.total);
A38b_bef(:,inp.Pg) = Pg38b_bef;
b38b = data.gen(:,10) * intmin;
[A38b,b38b] = time_relate(steps,A38b_bef,A38b,b38b);

A38 = [A38a; A38b];
b38 = [b38a; b38b];

ineq(38).A = A38;
ineq(38).b = b38;

%% Constraint 39
Xl39 = ones(1,len.Xl) .* data.load(:,4)' .* data.load(:,7)';
Xl39_bef = -Xl39;
Xg39 = -ones(1,len.Xg) .* data.gen(:,6)' .* data.gamma;
Xessd39 = -ones(1,len.Xessd) .* data.ess(:,15)' .* data.gamma;

A39 = zeros(1,len.total);
A39(:,inp.Xl) = Xl39;
A39(:,inp.Xg) = Xg39;
A39(:,inp.Xessd) = Xessd39;

A39_bef = zeros(1,len.total);
A39_bef(:,inp.Xl) = Xl39_bef;

b39 = zeros(size(A39,1),1);

[A39,b39] = time_relate(steps,A39_bef,A39,b39); 
ineq(39).A = A39;
ineq(39).b = b39;

%% CONNECTIVITY AND SEQUENCING CONSTRAINTS
%% Constraint 40
Xg40 = eye(len.Xg);
Sn40 = zeros(len.Xg,len.Sn);
for i = 1:len.Xg % Fills corresponding node that contains generator
    Sn40(i,data.gen(i,2)) = -1;
end

A40 = zeros(len.Xg,len.total);
A40(:,inp.Xg) = Xg40;
A40(:,inp.Sn) = Sn40;
A40 = A40(data.statgen==2,:); % Only use the non black-start DG

b40 = zeros(size(A40,1),1);

ineq(40).A = concA(steps,A40);
ineq(40).b = concB(steps,b40);

%% Constraint 41
Xbr41 = eye(len.Xbr);
Sn41 = zeros(len.Xbr,len.Sn);
for i = 1:len.Xbr
    Sn41(i,data.branch(i,2)) = -1;
end
A41 = zeros(len.Xbr,len.total);
A41(:,inp.Xbr) = Xbr41;
A41(:,inp.Sn) = Sn41;
A41 = A41(data.statbr == 2,:); % Only use the switchable line
b41 = zeros(size(A41,1),1);

ineq(41).A = concA(steps,A41);
ineq(41).b = concB(steps,b41);

%% Constraint 42
Xbr42 = eye(len.Xbr);
Sn42 = zeros(len.Xbr,len.Sn);
for i = 1:len.Xbr
    Sn42(i,data.branch(i,3)) = -1;
end
A42 = zeros(len.Xbr,len.total);
A42(:,inp.Xbr) = Xbr42;
A42(:,inp.Sn) = Sn42;
A42 = A42(data.statbr == 2,:); % Only use the switchable line
b42 = zeros(size(A42,1),1);

ineq(42).A = concA(steps,A42);
ineq(42).b = concB(steps,b42);

%% Constraint 43
% C.43 -- A
Xbr43a = eye(len.Xbr);
Sn43a = zeros(len.Xbr,len.Sn);
for i = 1:len.Xbr
    Sn43a(i,data.branch(i,2)) = -1;
end
Aeq43a = zeros(len.Xbr,len.total);
Aeq43a(:,inp.Xbr) = Xbr43a;
Aeq43a(:,inp.Sn) = Sn43a;
Aeq43a = Aeq43a(data.statbr == 1,:); % Only use the non-switchable line
beq43a = zeros(size(Aeq43a,1),1);

% C.43 -- B
Xbr43b = eye(len.Xbr);
Sn43b = zeros(len.Xbr,len.Sn);
for i = 1:len.Xbr
    Sn43b(i,data.branch(i,3)) = -1;
end
Aeq43b = zeros(len.Xbr,len.total);
Aeq43b(:,inp.Xbr) = Xbr43b;
Aeq43b(:,inp.Sn) = Sn43b;
Aeq43b = Aeq43b(data.statbr == 1,:); % Only use the non-switchable line
beq43b = zeros(size(Aeq43b,1),1);

Aeq43 = [Aeq43a; Aeq43b];
beq43 = [beq43a; beq43b];

equ(43).Aeq = concA(steps,Aeq43);
equ(43).beq = concB(steps,beq43);

%% Constraint 44
Sn44 = zeros(len.Xg,len.Sn);
for i = 1:len.Xg
    Sn44(i,data.gen(i,2)) = 1;
end
Xg44 = -eye(len.Xg);
Aeq44 = zeros(len.Xg,len.total);
Aeq44(:,inp.Sn) = Sn44;
Aeq44(:,inp.Xg) = Xg44;
Aeq44 = Aeq44(data.statgen == 1,:); % Only use the black-start DG
beq44 = zeros(size(Aeq44,1),1);

equ(44).Aeq = concA(steps,Aeq44);
equ(44).beq = concB(steps,beq44);

%% Constraint 45
Xl45 = eye(len.Xl);
Sn45 = zeros(len.Xl,len.Sn);
for i = 1:len.Xl
    Sn45(i,data.load(i,2)) = -1;
end
A45 = zeros(len.Xl,len.total);
A45(:,inp.Xl) = Xl45;
A45(:,inp.Sn) = Sn45;
A45 = A45(data.statload == 2,:); % Only use the switchable load
b45 = zeros(size(A45,1),1);

ineq(45).A = concA(steps,A45);
ineq(45).b = concB(steps,b45);

%% Constraint 46
Xl46 = eye(len.Xl);
Sn46 = zeros(len.Xl,len.Sn);
for i = 1:len.Xl
    Sn46(i,data.load(i,2)) = -1;
end
Aeq46 = zeros(len.Xl,len.total);
Aeq46(:,inp.Xl) = Xl46;
Aeq46(:,inp.Sn) = Sn46;
Aeq46 = Aeq46(data.statload == 1,:); % Only use the non-switchable load
beq46 = zeros(size(Aeq46,1),1);

equ(46).Aeq = concA(steps,Aeq46);
equ(46).beq = concB(steps,beq46);

%% Constraint 47
Xbr47 = -eye(len.Xbr);
Xbr47_bef = eye(len.Xbr);

A47 = zeros(len.Xbr,len.total);
A47(:,inp.Xbr) = Xbr47;
A47 = A47(data.statbr == 2,:); % Only use the switchable line

A47_bef = zeros(len.Xbr,len.total);
A47_bef(:,inp.Xbr) = Xbr47_bef;
A47_bef = A47_bef(data.statbr == 2,:); % Only use the switchable line

b47 = zeros(size(A47,1),1);
[ineq(47).A,ineq(47).b] = time_relate(steps,A47_bef,A47,b47);

%% Constraint 48
Xg48 = -eye(len.Xg);
Xg48_bef = eye(len.Xg);

A48 = zeros(len.Xg,len.total);
A48(:,inp.Xg) = Xg48;

A48_bef = zeros(len.Xg,len.total);
A48_bef(:,inp.Xg) = Xg48_bef;

b48 = zeros(size(A48,1),1);
[ineq(48).A,ineq(48).b] = time_relate(steps,A48_bef,A48,b48);

%% Constraint 49
Xl49 = -eye(len.Xl);
Xl49_bef = eye(len.Xl);
A49 = zeros(len.Xl,len.total);
A49(:,inp.Xl) = Xl49;
A49 = A49(data.statload == 2,:); % Only use the switchable load

A49_bef = zeros(len.Xl,len.total);
A49_bef(:,inp.Xl) = Xl49_bef;
A49_bef = A49_bef(data.statload == 2,:); % Only use the switchable load

b49 = zeros(size(A49,1),1);
[ineq(49).A,ineq(49).b] = time_relate(steps,A49_bef,A49,b49);

%% Running the MILP
RunMILP;
toc;
AllPrint;
AllPlot;
