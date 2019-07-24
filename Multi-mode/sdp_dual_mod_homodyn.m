function [Rho, Delta,fidel,distrace] = sdp_dual_mod_homodyn(rhot,projX_LARGE,medidas)

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

Obs=projX_LARGE;


Prob = (medidas);

Delta = sdpvar(length(projX_LARGE),1,'full','real');

F=[F,Delta>=0];


for i=1:length(projX_LARGE)
    F=[F,trace(Rho*Obs{i})<=Prob(i)+Delta(i)];
    F=[F,trace(Rho*Obs{i})>=Prob(i)-Delta(i)];
    measure = trace(Rho*Obs{i});
%     F = [F,measure>=Prob(i)*(1-Delta(i))];
%     F = [F,measure<=Prob(i)*(1-Delta(i))];
end



% cost function
% E = sdpvar(1,1,'Real');
% F = [F,E>=0];
E = sum(sum(Delta));

ops = sdpsettings('solver','mosek','verbose',1);
% ops.mosek.MSK_IPAR_NUM_THREADS=6;
SOLUTION=optimize(F,E,ops);


disp('DEBUGGING');
problema = double(SOLUTION.problem);
disp(yalmiperror(problema));

Rho = value(Rho);
Delta = value(Delta);

fidel=fidelity(rhot,Rho);
distrace=dist_trace(rhot,Rho);

