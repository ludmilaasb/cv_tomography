function projX_LARGE = make_projX_theta_homodyn_shifted_dual(L,Q,delta_theta,fock_size,arrede)


projX_theta_new = make_projX_theta_homodyn_shifted(L,Q,delta_theta,fock_size,arrede);

    
    projXtheta_col = reshape(transpose(projX_theta_new),[],1);
%     
%     projX_LARGE = cell(length(projXtheta_col));
%     
%     for i = 1:length(projXtheta_col)
%         projX_LARGE{i} = kron(projX_theta_new{i},projX_theta_new{i});
%     end
    
    
    projX_LARGE = cell(length(projXtheta_col),length(projXtheta_col));

    for i = 1:length(projXtheta_col)
        for j =1:length(projXtheta_col)
            projX_LARGE{i,j} = kron(projX_theta_new{i},projX_theta_new{i});
        end
    end
    
end