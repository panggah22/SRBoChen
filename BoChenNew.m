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

%% Running the MILP
RunMILP;
toc;