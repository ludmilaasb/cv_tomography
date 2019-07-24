function phi_beta = disp_wavelet_widt(n,beta,x,q0,p0)



% author: Ludmila Botelho
% date: Febreary, 2019

% Imputs: number state n, position in space x, width beta, displacement q0 and p0 

% Create wave fucntions for Fock States

phi_beta = (beta.^(1/2)./((2.^n).*factorial(n).*(pi.^(1./2))).^(1/2)).*hermiteH(n,beta.*(x-q0)).*exp(-(beta.^2.*(x-q0).^2)/2+1i*p0*q0 -((1i.*p0.*q0)/2));
end
