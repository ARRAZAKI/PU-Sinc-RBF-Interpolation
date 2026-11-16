function z = f1(x,y)
% Exemple 1
% f1(x,y) = 1/2 * y * (cos(4*(x^2 + y - 1)))^4

z = 0.5 * y .* (cos(4*(x.^2 + y - 1))).^4;
end
