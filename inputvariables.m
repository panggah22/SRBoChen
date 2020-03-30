function inp = inputvariables(len)

% Input ranges

inp.Pl = 0 + (1:len.Pl);
inp.Pbr = inp.Pl(end) + (1:len.Pbr);
inp.Pg = inp.Pbr(end) + (1:len.Pg);
inp.Pessc = inp.Pg(end) + (1:len.Pessc);
inp.Pessd = inp.Pessc(end) + (1:len.Pessd);
inp.Ql = inp.Pessd(end) + (1:len.Ql);
inp.Qbr = inp.Ql(end) + (1:len.Qbr);
inp.Qg = inp.Qbr(end) + (1:len.Qg);
inp.Qessc = inp.Qg(end) + (1:len.Qessc);
inp.Qessd = inp.Qessc(end) + (1:len.Qessd);
inp.U = inp.Qessd(end) + (1:len.U);
inp.Eess = inp.U(end) + (1:len.Eess);
inp.Sn = inp.Eess(end) + (1:len.Sn);
inp.Xl = inp.Sn(end) + (1:len.Xl);
inp.Xbr = inp.Xl(end) + (1:len.Xbr);
inp.Xg = inp.Xbr(end) + (1:len.Xg);
inp.Xessc = inp.Xg(end) + (1:len.Xessc);
inp.Xessd = inp.Xessc(end) + (1:len.Xessd);

inp.cont = 1:inp.Eess(end);
inp.disc = inp.Sn(1):inp.Xessd(end);

end