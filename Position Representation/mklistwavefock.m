function psis = mklistwavefock(n)


% author: Ludmila Botelho
%  date: March, 2018

% This function makes a list of psi_n(x)

% n = Number of elements


syms x %This is to generate the symbolic function for Fock State



psin = wavfunFock(n,x);     % psin = (1./((2.^n).*factorial(n).*(pi^(1/2))).^(1/2)).*hermiteH(n,x).*exp(-(x.^2)/2);


% The list of wave functions
psis = cell(length(psin),1);


%Converting to function handle, so you can use integral!

for i=1:length(psin)
    psis{i} = matlabFunction(psin(i));
end
