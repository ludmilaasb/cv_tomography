function psi_alpha = coherent_Xrep(q0,p0,x)



% author: Ludmila Botelho
% date: January, 2019

% Imputs: coherent state alpha = 1/(2^1/2)*(q0+ip0), position in space x

% Create wave fucntions for Coherent States

psi_alpha = (pi.^(-1./4)).*exp(-(((x-q0).^2)/2) + 1i.*p0.*x -((1i.*p0.*q0)/2));
end
