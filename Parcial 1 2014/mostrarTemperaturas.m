function mostrarTemperaturas(T)
    fprintf('\nTemperaturas nodales:\n');
    for i=1:length(T)
        fprintf('  T%d = %.4f\n',i,T(i));
    end
end
