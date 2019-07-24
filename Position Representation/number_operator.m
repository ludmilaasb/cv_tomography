
% Size

% -----------------------------------

L = 100;       

Q = 5; P = 5;

qmin = -Q; qmax = Q; deltaq = (qmax - qmin)/L;
q = qmin:deltaq:(qmax-deltaq);

pmin = -P; pmax = P; deltap = (pmax - qmin)/L;
p = pmin:deltap:(pmax-deltap);
% --------------------------------------




% Make a list of psi_n(x)

% ----------------------------------------
n = 0:3; % Number of elements


syms x %This is to generate the symbolic function for Fock State

% psin = (1./((2.^n).*factorial(n).*(pi^(1/2))).^(1/2)).*hermiteH(n,x).*exp(-(x.^2)/2);

psin = wavfunFock(n,x);


% The list
psis = cell(length(psin),1);


%Converting to function handle

for i=1:length(psin)
    psis{i} = matlabFunction(psin(i));
end


% wnumproj = cell(length(n):1);
N  = zeros(L,L);
for i = 1:length(n)
%     wnumproj{i} = wignerf(psis{i},L,Q,P).*n(i);
    N = N + wignerf(psis{i},L,Q,P).*n(i);
end

% save(sprintf('N_%d_%d.mat',0,3),'N')

