function [F] = genF(xnode,Q)
    if (size(xnode,1) == 3)
        F = (Q*area_ele(xnode)/3)*ones(3,1);
    else
        F = (Q*area_ele(xnode)/4)*ones(4,1);
    end
end