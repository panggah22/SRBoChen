function [a,b] = time_relate(steps,Aeq_before,Aeq_after,beq)

con = cell(steps-1,steps);
conb = cell(steps-1,1);
for i = 1:steps-1
    con{i,i} = Aeq_before;
    con{i,i+1} = Aeq_after;
    conb{i} = beq;
end
for i = 1:steps-1
    for j = 1:steps
        if isempty(con{i,j})
            con{i,j} = zeros(size(Aeq_before));
        end
    end
end
a = cell2mat(con);
b = cell2mat(conb);

end