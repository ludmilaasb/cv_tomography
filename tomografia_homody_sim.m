% Autor: Ludmila Botelho

% ludmilaasb@gmail.com 

% I've made this function to simulate a homodyne tomography 
% and reconstructy a density matrix RHO.
% The imputs are:  L (how many points we are going to measure), 
% Q (the lenght of the space of measurement)
% fock_size (how many elements of fock space are we going to use)
% and the density_matriz we are going to evaluate
% The outputs are: Rho (the reconstructed state)
% Delta (the error, it's interesting if we add noise
% to the density_matriz)
% distance_trace (the trace distance between Rho and 
% density_matrix)
% and fidel (the fidelity F(density_matrix,Rho)


% It is important to note that you should have
% a mosek license (or any other lincense to make
% the yalmip works)

% U can ALWAYS send me and email! I'm not
% a programmer, I just did it to help to solve
% my stuff!

% Enjoy!


function [Rho,Delta,distance_trace,fidel]=tomografia_homody_sim(L,Q,delta_theta,fock_size,density_matrix)

% L = 100;       %State Size
% Q = 5;

%              
%------------------------------------------------
qmin = -Q; qmax = Q; deltaq = (qmax - qmin)/L;
q = qmin:deltaq:(qmax-deltaq);
%---------------------------------------------------

theta = deg2rad(0:delta_theta:180);

n = fock_size;
ns = 0:fock_size;

% wavfunFock(ns(1),-5);


% base_c = basecanon(n);

a = annihilfock_mobral(n);
N = a'*a;



% ~ coherent ~

rhot = density_matrix;



% syms x
% 
% psisn = wavfunFock(ns,-4);

% a = annihilfock_mobral(n);
% 
% projX = cell(length(q),1);
% 
braX = zeros(length(q),length(ns));

for j=1:length(braX)
    for k=1:length(ns)
        
        braX(j,k)= wavfunFock(ns(k),q(j));
        
    end
end


U = cell(length(theta),1);

projX_theta = cell(length(q),length(theta));

medidas = zeros(size(projX_theta));


for k =1:length(q)
    for l = 1:length(theta)
        
        U{l} = zeros(size(N));
        U{l} = unit_rot_theta(N,theta(l));

        projX_theta{k,l} = zeros(length(ns),length(ns));
        projX_theta{k,l} = U{l}'*braX(k,:)'*braX(k,:)*U{l};
        
        medidas(k,l) = trace(rhot*projX_theta{k,l});
             
%     projX{k} = zeros(length(ns),length(ns));
%     projX{k} = ketX(k,:)'*ketX(k,:);

    end
end

noisemin=-0.2;
noisemax=0.2;
noise=noisemin+rand(length(q),length(theta))*(noisemax-noisemin);




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

Obs=projX_theta;


Prob = (medidas+noise);

Delta = sdpvar(length(q),length(theta),'full','real');

F=[F,Delta>=0];


for i=1:length(q)
    for j=1:length(theta)
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

distance_trace=dist_trace(rhot,Rho);
fidel = fidelity(rhot,Rho);

end
