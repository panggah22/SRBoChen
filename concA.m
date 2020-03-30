function fuu = concA(SetSteps,in)

ce = cell(SetSteps);
for i = 1:SetSteps
    for j = 1:SetSteps
        if i == j
            ce{i,j} = in;
        else
            ce{i,j} = zeros(size(in));
        end
    end
end
fuu = cell2mat(ce);

end