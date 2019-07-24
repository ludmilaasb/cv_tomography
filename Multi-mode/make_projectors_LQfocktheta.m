function projX_theta = make_projectors_LQfocktheta(L,Q,n,delta_theta)


% L = Intervalo do espaço
% Q = Extensão
% n = número do espaço de Fock
% 

%              _Intervalos_
%------------------------------------------------
qmin = -Q; qmax = Q; deltaq = (qmax - qmin)/L;
q = qmin:deltaq:(qmax-deltaq);
%---------------------------------------------------

theta = deg2rad(0:delta_theta:180);


ns = 0:n;



% ----------- Ket X ------------------------

ketX = zeros(length(ns),length(q));

mketx2 = cell(length(q),length(q));


for i=1:length(ns)
    for j=1:length(q)
        ketX(i,j)= wavfunFock(ns(i),q(j));
    end
end
    
           
       
         
for j=1:length(q)
    for k=1:length(q)
        mketx2{j,k} = zeros(length(ns)^2,1);
        mketx2{j,k} = reshape((ketX(:,j)*ketX(:,k)')',[],1);
    end
end
      
ketXbase2 = reshape(mketx2',[],1);


% ------------ Unitárias e^(i*theta*N)-------

mU = cell(length(theta),length(theta));


for i = 1:length(theta)
    for j =1:length(theta)
%         mU{i,j} = zeros(length(N)^2,length(N)^2);
        mU{i,j} = kron(unit_rot_theta(N,theta(i)),unit_rot_theta(N,theta(j)));
    end
end

U = reshape(mU',[],1);


% ------------ Projetores-------------------

projX_theta = cell(length(ketXbase2),length(U));

for i=1:length(ketXbase2)
    for j=1:length(U)
        projX_theta{i,j} = U{j}'*ketXbase2{i}*ketXbase2{i}'*U{j};
    end
end

