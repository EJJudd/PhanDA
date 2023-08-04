function CE = coefefficiency(X1, X2)

CE = 1 - sum((X1 - X2).^2) / sum((X1 - mean(X1)).^2);

end

