function phi_beta_n = eign_harmonic(n,beta,x)



% author: Ludmila Botelho
% date: Febreary, 2019

% Imputs: number state n, position in space x

% Create wave fucntions for Fock States

phi_beta_n = ((beta.^(1/2)./(2.^n).*factorial(n).*(pi.^(1./2))).^(1/2)).*hermiteH(n,beta.*x).*exp(-(beta.^2.*x.^2)/2);
end
