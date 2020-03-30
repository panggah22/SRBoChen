function fuu = lasts(steps,in)

con = cell(1,steps);
for i = 1:steps
    if i == steps
        con{i} = in;
    else 
        con{i} = zeros(size(in));
    end
end
fuu = cell2mat(con);

end