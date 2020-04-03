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

%% Print the MILP Result
% load energization
fprintf('-------------------------\n');
fprintf('  Energization of loads\n');
fprintf('-------------------------\n');
fprintf('Timestep | Load energized\n');
fprintf('-------------------------\n');
for i = 1:steps
    fprintf('  %d         ',i);
    for j = 1:size(ts_load{i,2},1)
        if ~isempty(ts_load{i,2})
            if j == size(ts_load{i,2},1)
                fprintf('L%d',ts_load{i,2}(j));
            else
                fprintf('L%d, ',ts_load{i,2}(j));
            end
        else
            fprintf('     ');
        end
    end
    fprintf('\n');
end

% DG energization
fprintf('--------------------------\n');
fprintf('Energization of generators\n');
fprintf('--------------------------\n');
fprintf('Timestep | DG energized\n');
fprintf('--------------------------\n');
for i = 1:steps
    fprintf('  %d         ',i);
    for k = 1:size(ts_gen{i,2},1)
        if ~isempty(ts_gen{i,2})
            if k == size(ts_gen{i,2},1)
                fprintf('DG%d',ts_gen{i,2}(k));
            else
                fprintf('DG%d, ',ts_gen{i,2}(k));
            end
        else
            fprintf('     ');
        end
    end
    fprintf('\n');
end

% load energization
fprintf('---------------------------\n');
fprintf('  Energization of branches\n');
fprintf('---------------------------\n');
fprintf('Timestep | Branch energized\n');
fprintf('---------------------------\n');
for i = 1:steps
    fprintf('  %d         ',i);
    for j = 1:size(ts_branch{i,2},1)
        if ~isempty(ts_branch{i,2})
            if j == size(ts_branch{i,2},1)
                fprintf('B%d',ts_branch{i,2}(j));
            else
                fprintf('B%d, ',ts_branch{i,2}(j));
            end
        else
            fprintf('     ');
        end
    end
    fprintf('\n');
end

% % load energization
% fprintf('--------------------------\n');
% fprintf('   Energization of ESS\n');
% fprintf('--------------------------\n');
% fprintf('Timestep | Nodes energized\n');
% fprintf('--------------------------\n');
% for i = 1:steps
%     fprintf('  %d         ',i);
%     for j = 1:size(ts_node{i,2},1)
%         if ~isempty(ts_node{i,2})
%             if j == size(ts_node{i,2},1)
%                 fprintf('ESS%d',ts_node{i,2}(j));
%             else
%                 fprintf('ESS%d, ',ts_node{i,2}(j));
%             end
%         else
%             fprintf('     ');
%         end
%     end
%     fprintf('\n');
% end