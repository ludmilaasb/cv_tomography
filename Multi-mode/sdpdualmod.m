clear;


load('HermG1Fock21x21.mat');


n =13;
ns = 0:n;

truncamento = (length(ns))^2;

rhot = zeros(truncamento,truncamento);

for i=1:(truncamento)
    for j=1:(truncamento)
        rhot(i,j) = HermGFock(i,j);
    end
end



base_c = basecanon(n);




mbase2 = cell(length(ns),length(ns));

for i=1:length(ns)
    for j=1:length(ns)
        mbase2{i,j} = zeros(length(ns)^2,1);
        mbase2{i,j} = reshape((base_c{i}*base_c{j}')',[],1);
        base2 = reshape(mbase2',[],1);

    end
end

projectors = cell(length(ns)^2,length(ns)^2);

for i = 1:length(ns)^2
    for j = 1:length(ns)^2
        projectors{i,j} = base2{i}*base2{j}';
    end
end

medidas =  zeros(size(rhot));

for i = 1:length(rhot)
    for j = 1:length(rhot)
%         medidas(i,j) = base2{i}'*rhot*base2{j};
        medidas(i,j) = trace(rhot*projectors{i,j});
    end
end

x = zeros(size(rhot)); sigma=.3;mu=0;
for i=1:length(rhot)
    for j=1:length(rhot)
        z0 = (sqrt(-2.0 * log(rand)) * cos(2*pi * rand));
        x(i,j) = sigma*z0 +mu;
    end
end
histogram(x)


noise = x;

df = length(rhot);
% Nobs = numel(projectors);
% cleaning yalmip memory
yalmip('clear');

F = class('double');

% defining the SDP variables
Rho = sdpvar(df,df,'hermitian','complex');

% standard constraints
F=[Rho>=0];
F=[F,trace(Rho)==1];

% observables

Obs=projectors;


Prob = (HermGFock +noise);

Delta = sdpvar(df,df,'full','real');

F=[F,Delta>=0];


for i=1:df
    for j=1:df
        F=[F,trace(Rho*Obs{i,j})<=Prob(i,j)+Delta(i,j)];
        F=[F,trace(Rho*Obs{i,j})>=Prob(i,j)-Delta(i,j)];
        measure = trace(Rho*Obs{i,j});
%     F = [F,measure>=Prob(i)*(1-Delta(i))];
%     F = [F,measure<=Prob(i)*(1-Delta(i))];
    end
end



% cost function
% E = sdpvar(1,1,'Real');
% F = [F,E>=0];
E = sum(sum(Delta));

ops = sdpsettings('solver','mosek','verbose',1);
ops.mosek.MSK_IPAR_NUM_THREADS=6;
SOLUTION=optimize(F,E,ops);


disp('DEBUGGING');
problema = double(SOLUTION.problem);
disp(yalmiperror(problema));

Rho = value(Rho);
Delta = value(Delta);

fidelity(rhot,Rho)
