clear

L = 100;

n = 0:30;
Q = 5;
P = 5;

% -----------------------------------------------
qmin = -Q; qmax = Q; deltaq = (qmax - qmin)/L;
q = qmin:deltaq:(qmax-deltaq);

pmin = -P; pmax = P; deltap = (pmax - qmin)/L;
p = pmin:deltap:(pmax-deltap);
% -----------------------------------------------

alpha1 = 1;
alpha2 = -1;

load('wprojectors_0_30.mat'); 



walpha1 = zeros(L,L);
walpha2 = zeros(L,L);


for i = 1:length(n)
    for j = 1:length(n)
        walpha1 = walpha1 + exp(-abs(alpha1).^2)*((alpha1.^(n(i))*alpha1^(n(j)))/(((factorial(n(i)))^1/2)*((factorial(n(j)))^1/2)))*(wprojectors{i,j});
        walpha2 = walpha2 + exp(-abs(alpha2).^2)*((alpha2.^(n(i))*alpha2^(n(j)))/(((factorial(n(i)))^1/2)*((factorial(n(j)))^1/2)))*(wprojectors{i,j});
    end
end





figure(1)
s1 = surf(p,q,walpha1);
xlabel('Q')
ylabel('P')
zlabel('Feiura')
s1.EdgeColor = 'none';

figure(2)
s2= surf(p,q,walpha2);
xlabel('Q')
ylabel('P')
zlabel('Feiura')
s2.EdgeColor = 'none';
