function [eG] = assembly(eG,eL,ele)
    if (size(eG,2) == 1)
        eG(ele) = eG(ele) + eL;
    else
        eG(ele,ele) = eG(ele,ele) + eL;
    end
end