function len = lengthbochen(data)

% Continuous variables
len.Pl = data.num_load; len.Pbr = data.num_branch; len.Pg = data.num_gen;
len.Pessc = data.num_ess; len.Pessd = len.Pessc;
len.Ql = len.Pl; len.Qbr = len.Pbr; len.Qg = len.Pg;
len.Qessc = len.Pessc; len.Qessd = len.Pessd; len.U = data.num_bus;
len.Eess = data.num_ess;

% Discrete and binary variables
len.Sn = data.num_bus; len.Xl = len.Pl; len.Xbr = data.num_branch;
len.Xg = data.num_gen; len.Xessc = data.num_ess; len.Xessd = len.Xessc;

len.total = len.Pl + len.Pbr + len.Pg + len.Pessc + len.Pessd ...
    + len.Ql + len.Qbr + len.Qg + len.Qessc + len.Qessd + len.U + len.Eess...
    + len.Sn + len.Xl + len.Xbr + len.Xg + len.Xessc + len.Xessd;
len.Ptot = len.Pl + len.Pbr + len.Pg + len.Pessc + len.Pessd;
len.Qtot = len.Ql + len.Qbr + len.Qg + len.Qessc + len.Qessd;
len.U = len.U;
len.Eess = len.Eess;
len.discrete = len.Sn + len.Xl + len.Xbr + len.Xg + len.Xessc + len.Xessd;

end