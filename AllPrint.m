Timestep = 1:steps;
xnode = x(inp.Sn,:); 
xload = x(inp.Xl,:);
xgen = x(inp.Xg,:);
xbranch = x(inp.Xbr,:);

ts_node = cell(steps,2);
ts_load = cell(steps,2);
ts_gen = cell(steps,2);
ts_branch = cell(steps,2);


for i = 1:steps
    ts_node{i,1} = data.bus(xnode(:,i) == 1,2);
    ts_load{i,1} = data.load(xload(:,i) == 1,3);
    ts_gen{i,1} = data.gen(xgen(:,i) == 1,1);
    ts_branch{i,1} = data.branch(xbranch(:,i) == 1,1);
    if i > 1
        ts_node{i,2} = setdiff(ts_node{i,1},ts_node{i-1,1});
        ts_load{i,2} = setdiff(ts_load{i,1},ts_load{i-1,1});
        ts_gen{i,2} = setdiff(ts_gen{i,1},ts_gen{i-1,1});
        ts_branch{i,2} = setdiff(ts_branch{i,1},ts_branch{i-1,1});
    else
        ts_node{i,2} = ts_node{i,1};
        ts_load{i,2} = ts_load{i,1};
        ts_gen{i,2} = ts_gen{i,1};
        ts_branch{i,2} = ts_branch{i,1};
    end
end