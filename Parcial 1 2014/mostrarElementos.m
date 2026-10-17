function mostrarElementos(Ke,Fe)
    for e = 1:length(Ke)
        fprintf('\n-- Elemento e%d --\n',e);
        fprintf('Matriz de rigidez:\n'); disp(Ke{e});
        fprintf('Vector de cargas:\n'); disp(Fe{e});
    end
end
