function wprojectors = mklistfckproj(n,L,Q,P)

% author: Ludmila Botelho
% date: March, 2018


% This function makes wigner functions of fock projectors.

% IMPUTS:

% n = number of fock basis
% L = resolution
% Q,P = limits of the quadrature space 


% Creates the a list of wave functions so we can Wigner it: 
psis = mklistwavefock(n);



% Size

% -----------------------------------

qmin = -Q; qmax = Q; deltaq = (qmax - qmin)/L;
q = qmin:deltaq:(qmax-deltaq);

pmin = -P; pmax = P; deltap = (pmax - qmin)/L;
p = pmin:deltap:(pmax-deltap);
% --------------------------------------



wprojectors = cell(length(n):length(n));

for i = 1:length(n)
    for j = 1:length(n)
        wprojectors{i,j} = zeros(L:L);
       func = @(x,q,p) (1/(2*pi))*psis{i}(q-x/2).*conj(psis{j}(q+x/2)).*exp(1i*p*x);
       
       for k = 1:L
           for l = 1:L
                wprojectors{i,j}(k,l) = real(integral(@(x)func(x,q(k),p(l)),-Inf,Inf));
               
           end
       end

    end
end
