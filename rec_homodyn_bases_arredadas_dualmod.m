function [arrede,rec_rhos] = rec_homodyn_bases_arredadas_dualmod(L,Q,delta_theta,Lim,fock_size,Arr_int,medidas)



arrede=0:Arr_int:Lim;

rec_rhos = cell(length(arrede),1);


% [xx,yy] = size(medidas);
[xx,yy] = size(medidas);
for ii =1:length(arrede)
    
    projX_LARGE = make_projX_theta_homodyn_shifted_dual(L,Q,delta_theta,fock_size,arrede(ii));
    
    
    
    
    
%     projX_theta_new = make_projX_theta_homodyn_shifted(L,Q,delta_theta,fock_size,arrede(ii));
% 
%     
%     projXtheta_col = reshape(transpose(projX_theta_new),[],1);
%     
%     projX_LARGE = cell(length(projXtheta_col),length(projXtheta_col));
% 
%     for i = 1:length(projXtheta_col)
%         for j =1:length(projXtheta_col)
%             projX_LARGE{i,j} = kron(projX_theta_new{i},projX_theta_new{j});
%         end
%     end
    
    rec_rhos{ii} = zeros(fock_size^2,fock_size^2);
    
      for i=1:xx
          rec_rhos{ii} =rec_rhos{ii} + medidas(i)*projX_LARGE{i};
      end
      
      
%     for i=1:xx
%         for j=1:yy
%             rec_rhos{ii} =rec_rhos{ii} + medidas(i,j)*projX_LARGE{i,j};
%         end
%     end
    
end

end

