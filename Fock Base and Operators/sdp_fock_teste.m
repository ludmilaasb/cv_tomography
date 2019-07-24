

% cleaning yalmip memory
yalmip('clear');

F = class('double');

% defining the SDP variables
Rho = sdpvar(df, df, 'hermitian', 'complex');

% standard constraints
F = [F,Rho>=0];
F = [F,trace(Rho)==1];

% observables

Obs=cell(Nobs,1);


Prob = size(Nobs);

Delta = sdpvar(Nobs,'Real');

F = [F,Delta(:)>=0];

for i=1:Nobs
    F = [F,trace(Rho*Obs{i})<=Prob(i)+Delta(i)];
    F = [F,trace(Rho*Obs{i})>=Prob(i)-Delta(i)];
%     measure = trace(Rho*Obs{i});
%     F = [F,measure>=Prob(i)*(1-Delta(i))];
%     F = [F,measure<=Prob(i)*(1-Delta(i))];
end



% cost function
% E = sdpvar(1,1,'Real');
% F = [F,E>=0];
E = sum(Delta);



SOLUTION = optimize(F,E);
disp('DEBUGGING');
problema = double(SOLUTION.problem);
disp(yalmiperror(problema));

Rho = value(Rho);
Delta = value(Delta);