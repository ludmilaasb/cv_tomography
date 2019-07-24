
% Make a list of psi_n(x)

n = 0:1; % Number of elements


syms x %This is to generate the symbolic function for Fock State

% psin = (1./((2.^n).*factorial(n).*(pi^(1/2))).^(1/2)).*hermiteH(n,x).*exp(-(x.^2)/2);

psin = wavfunFock(n,x);


% The list
psis = cell(length(psin),1);


%Converting to function handle

for i=1:length(psin)
    psis{i} = matlabFunction(psin(i));
end


% Size

% -----------------------------------

L = 100;       

Q = 5; P = 5;

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
