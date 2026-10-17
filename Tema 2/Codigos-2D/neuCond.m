function [F] = neuCond(F,cN,xnode)
    for i=1:size(cN,1)
        x=xnode(cN(i,1),:)-xnode(cN(i,2),:); % Side Coordinaates
        l = sqrt(x*transpose(x));                      % Side Size
        lado = [cN(i,1) cN(i,2)];            % Side Nodes
        q = cN(i,3);
        F(lado) = F(lado) - (q*l/2)*ones(2,1);     
    end
end