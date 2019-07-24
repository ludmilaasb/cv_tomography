clear;

rhot = RandomDensityMatrix(21);

base_c=basecanon(20);

projectors = cell(length(base_c):length(base_c));

for i=1:length(base_c)
    for j=1:length(base_c)
%         projectors{i,j}=zeros(size(rhot));
        projectors{i,j}=base_c{i}*base_c{j}';
    end
end

medidas = zeros(length(rhot),length(rhot));

for i=1:length(rhot)
    for j=1:length(rhot)
        medidas(i,j) = trace(rhot*projectors{i,j});
    end
end

df = length(rhot);
Nobs = numel(projectors);
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


Prob = medidas;

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