function fuu = concB(SetSteps,in)

ce = cell(SetSteps,1);
for i = 1:SetSteps
    ce{i} = in;
end
fuu = cell2mat(ce);

end