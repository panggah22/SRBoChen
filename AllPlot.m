%% Plot figures
%% ESS
figure
ESSplot = x(inp.Eess,:)./data.ess(:,4)'*100;
plot(ESSplot);
ylim([9 40])

%% Active power generated
esscdp = [x(inp.Pessc,:); x(inp.Pessd,:)]';
Pesscd = zeros(size(esscdp,1),1);
for i = 1:size(esscdp,1)
    if esscdp(i,1) == 0
        Pesscd(i) = esscdp(i,2);
    else
        if esscdp(i,2) == 0
            Pesscd(i) = -esscdp(i,1);
        else 
            Pesscd(i) = 0;
        end
    end
end
figure
P_all = [Pesscd x(inp.Pg,:)']*data.MVAbase*1000;

bar(P_all,'stacked');

%% Rective power generated
esscdq = [x(inp.Qessc,:); x(inp.Qessd,:)]';
Qesscd = zeros(size(esscdq,1),1);
for i = 1:size(esscdq,1)
    if esscdq(i,1) == 0
        Qesscd(i) = esscdq(i,2);
    else
        if esscdp(i,2) == 0
            Qesscd(i) = -esscdq(i,1);
        else 
            Qesscd(i) = 0;
        end
    end
end
figure
Q_all = [Qesscd x(inp.Qg,:)']*data.MVAbase*1000;

bar(Q_all,'stacked');
