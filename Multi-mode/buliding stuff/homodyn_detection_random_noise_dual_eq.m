function [medidas,projX_LARGE,noise] = homodyn_detection_random_noise_dual_eq(L,Q,delta_theta,fock_size,density_matrix,noisemin,noisemax)

% clear;
% 
% L=10;
% Q=5;
% delta_theta=60;
% 
% fock_size = 5;

% qmin = -Q; qmax = Q; deltaq = (qmax - qmin)/L;
% q = qmin:deltaq:(qmax-deltaq);
% 
% theta = deg2rad(0:delta_theta:180);

projX_theta1 = make_projX_theta_homodyn_shifted(L,Q,delta_theta,fock_size,0);


%     projXtheta_col = reshape(transpose(projX_theta1),[],1);
%     
%     projX_LARGE = cell(length(projXtheta_col),1);
%     
%     for i = 1:length(projXtheta_col)
%         projX_LARGE{i} = kron(projX_theta1{i},projX_theta1{i});
%     end
%     
%     
% medidas = size(projX_LARGE);
% 
% for i = 1:length(projXtheta_col)
%     medidas(i) = trace(density_matrix*projX_LARGE{i});
% end
% -------------------------------------------------------------------
projXtheta_col = reshape(transpose(projX_theta1),[],1);

projX_LARGE = cell(length(projXtheta_col),length(projXtheta_col));




for i = 1:length(projXtheta_col)
    for j =1:length(projXtheta_col)
        projX_LARGE{i,j} = kron(projX_theta1{i},projX_theta1{j});
    end
end



medidas = size(projX_LARGE);

for i = 1:length(projXtheta_col)
    for j = 1:length(projXtheta_col)
        medidas(i,j) = trace(density_matrix*projX_LARGE{i,j});
    end
end


noise=noisemin+rand(size(medidas))*(noisemax-noisemin);

medidas = medidas + noise;

end