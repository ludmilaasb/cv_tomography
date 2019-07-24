function psi_n = wavfunFock(n,x)



% author: Ludmila Botelho
% date: January, 2018

% Imputs: number state n, position in space x

% Create wave fucntions for Fock States

psi_n = (1./((2.^n).*factorial(n).*(pi.^(1./2))).^(1/2)).*hermiteH(n,x).*exp(-(x.^2)/2);
end
