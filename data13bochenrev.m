function data = data13bochen
% Adding initial comments
%% General data for conversion
data.kVbase = 12.6;
data.MVAbase = 100;
data.zbase = (data.kVbase^2)/data.MVAbase;

%% Load data
data.bus = [
%  Load_no  P(kW)   Q(kW)   Weight      S_U     S_D     D_L     alpha   Status
    632     100     40      1.62        2       1.2     2.2     0.7     2
    634     100     40      7.94        1.8     1.1     2.3     0.8     2
    645     200     100     3.11        2.4     1.1     1.1     1.0     2
    646     200     100     8.29        2.5     1.1     3.6     0.7     2
    671     300     180     1.66        2.3     1.0     2.7     1.3     2
    692     200     150     6.02        2.7     1.0     1.2     0.6     2
    675     200     40      2.63        3.4     1.3     2.7     0.7     1
    611     200     120     1.54        2.1     1.3     2.1     0.6     2
    652     40      30      3.89        3.0     1.1     2.0     0.9     0
    650     0       0       0           0       0       0       0       2
    633     0       0       0           0       0       0       0       2
    680     0       0       0           0       0       0       0       2
    684     0       0       0           0       0       0       0       2
    ];

% Converting loads to Per Unit value
data.bus(:,2:3) = data.bus(:,2:3)./(data.MVAbase*1000);
[data.num_bus,~] = size(data.bus);
data.bus = [(1:data.num_bus)' data.bus]; % for internal numbering

data.load = data.bus(data.bus(:,3)>0,:);
data.num_load = size(data.load,1);
data.load = [(1:data.num_load)' data.load];

%% Branch data
data.branch = [
%   From    To      Length  Capacity    Status
    650     632     2000    1500        2
    632     633     500     1000        2
    633     634     200     500         2
    632     645     500     1000        2
    645     646     300     800         2
    632     671     2000    1500        2
    671     692     10      800         0
    692     675     500     800         2
    671     684     300     800         0
    684     611     300     800         1
    684     652     800     800         1
    671     680     1000    1500        2
    633     692     2000    1000        2
    646     611     2000    800         2
    675     680     1500    1000        2
    ];

%Converting impedance to Per Unit value
milebranch = data.branch(:,3)./5280;
impedance = [milebranch*0.2     milebranch*0.8];
impedance = impedance./data.zbase;
data.branch = [data.branch(:,1:2)   impedance   data.branch(:,3:end)];
[data.num_branch,~] = size(data.branch);
fromto = zeros(data.num_branch,2);
for i = 1:data.num_branch % Internal numbering
    fromto(i,1) = data.bus((data.bus(:,2) == data.branch(i,1)),1);
    fromto(i,2) = data.bus((data.bus(:,2) == data.branch(i,2)),1);
end
data.branch = [(1:data.num_branch)' fromto data.branch]; % for internal numbering
data.branch(:,9) = data.branch(:,9)./(data.MVAbase*1000);

%% Generator data
% Status list
% 1 = Substation / Black-start
% 2 = Non-black-start
% 0 = in Fault
data.gen = [
%   Node    PF      FRR     Pg_max      Pg_min  Qg_max  Qg_min  Pramp   Status
    650     1       5       12000       0       9000    -8000   1000    1
    646     0.8     5       800         50      500     -400    200     2
    680     0.8     5       1000        100     1000    -800    300     2
    ];

% Converting generations to Per Unit value
data.gen(:,4:8) = data.gen(:,4:8)./(data.MVAbase*1000);
data.gen(:,3) = data.gen(:,3)./100; % Change to percentage
[data.num_gen,~] = size(data.gen);
n_gen = zeros(data.num_gen,1);
for i = 1:data.num_gen
    n_gen(i,1) = data.bus((data.bus(:,2) == data.gen(i,1)),1);
end
data.gen = [(1:data.num_gen)' n_gen data.gen];

%% ESS data
data.ess = [
%   Node    Rating  Socini  Socmin  Socmax  C_eff   D_eff   PCmin   PCmax   QCmin   QCmax   PDmin   PDmax   QDmin   QDmax   PRC     PRD     QRC     QRD
    632     200     30      10      100     0.9     0.9     0       5       0       2       0       5       0       2       250     250     250     250
    ];

% Converting power to Per Unit Value
data.ess(:,[2 (16:19)]) = data.ess(:,[2 (16:19)])./(data.MVAbase*1000);
data.ess(:,8:15) = data.ess(:,8:15)./data.MVAbase;
data.ess(:,3:5) = data.ess(:,3:5)./100; % Change into percentage
[data.num_ess,~] = size(data.ess);
n_ess = zeros(data.num_ess,1);
for i = 1:data.num_ess
    n_ess(i,1) = data.bus((data.bus(:,2) == data.ess(i,1)),1);
end
data.ess = [(1:data.num_ess) n_ess data.ess];

%% Other data
data.gamma = 0.05;
data.statgen = data.gen(:,end);
data.statload = data.load(:,end);
data.statbr = data.branch(:,end);
data.statbus = data.bus(:,end);
end