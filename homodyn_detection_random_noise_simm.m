function [medidas,projX_theta,noise] = homodyn_detection_random_noise_simm(L,Q,delta_theta,fock_size,density_matrix,noisemin,noisemax)

%              _Intervalos_
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

% noisemin=-0.2;
% noisemax=0.2;
noise=noisemin+rand(length(q),length(theta))*(noisemax-noisemin);

medidas = medidas + noise;
end