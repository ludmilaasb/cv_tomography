function [arrede,rec_rhos] = rec_homodyn_bases_arredadas(L,Q,delta_theta,D,Arr_int,medidas)


% author: Ludmila Botelho
% date: February, 2019

% Imputs: Dimension of space -> D
%         "Arrede"interval -> Arr_int
%         Measuring Results -> medidas

%  It is an aproximated reconstruction off the density matrix state through
%  a shifted fock basis. This version I will test homodyne tomography.




arrede=0:Arr_int:D;

rec_rhos = cell(length(arrede),1);


[qs,thetas] = size(medidas);
for ii =1:length(arrede)
    
    projX_theta_new = make_projX_theta_homodyn_shifted(L,Q,delta_theta,D,arrede(ii));

    rec_rhos{ii} = zeros(D,D);
    for i=1:qs
        for j=1:thetas
            rec_rhos{ii} =rec_rhos{ii} + medidas(i,j)*projX_theta_new{i,j};
        end
    end
    
    
    
    
end
    

    
        



    
      
                




