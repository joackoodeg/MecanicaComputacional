function [K,F] = robCond(K,F,cR,xnode)
    for i=1:size(cR,1)
        x=xnode(cR(i,1),:)-xnode(cR(i,2),:); % Side Coordinaates
        l = sqrt(x*transpose(x));                  % Side Size
        lado = [cR(i,1) cR(i,2)];            % Side Nodes
        hT = cR(i,3)*cR(i,4);
        K(lado,lado) = K(lado,lado)+cR(i,3)*l/6*[2 1;1 2];
        F(lado) = F(lado) + (l*hT/2)*ones(2,1);     
    end
end