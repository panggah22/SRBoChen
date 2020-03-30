function a = initials(SetSteps,Aeq)

conA = cell(1,SetSteps);
for i = 1:SetSteps
    if i == 1
        conA{i} = Aeq;
    else
        conA{i} = zeros(size(Aeq));      
    end
end
a = cell2mat(conA);

end