
function projX_theta_new = make_projX_theta_homodyn_shifted(L,Q,delta_theta,D,arrede)






n_arredado = arrede:(arrede+D-1);

qmin = -Q; qmax = Q; deltaq = (qmax - qmin)/L;
q = qmin:deltaq:(qmax-deltaq);

theta = deg2rad(0:delta_theta:180);



    projX_theta_new = cell(length(q),length(theta));
    for i =1:length(q)
        for j = 1:length(theta)
            projX_theta_new{i,j} = zeros(D,D);
        end
    end
    

    psis = zeros(length(q),length(n_arredado));
    
    for j=1:length(q)
        for k=1:length(n_arredado)
            psis(j,k)= wavfunFock(n_arredado(k),q(j));
        end
    end
    
    
    
    for i =1:length(q)
        for j = 1:length(theta)
            for k = 1:length(n_arredado)
                for l = 1:length(n_arredado)
                    projX_theta_new{i,j}(k,l)= conj(psis(i,k))*(psis(i,l))*exp(1i*(n_arredado(k)-n_arredado(l))*theta(j));
                end
            end
        end
    end
    
    
    
end

    