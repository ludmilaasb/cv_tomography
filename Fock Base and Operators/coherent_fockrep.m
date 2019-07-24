function ket_alpha = coherent_fockrep(alpha,d)

% d = Number of photons spaces dimension
% alpha = 1/sqrt(2) *(q + ip)^2

ket_alpha = zeros(d+1,1);



for n=1:(length(ket_alpha))
    ket_alpha(n)= exp(-(1/2)*abs(alpha)^2)*(alpha^(n-1)/sqrt(factorial(n-1)));
end

    
    
%     ket_alpha = exp(-(1/2)*abs(alpha)^2)*base_c{1};
%     ket_alpha = ket_alpha + exp(-(1/2)*abs(alpha)^2)*(alpha^(n)/sqrt(factorial(n)))*(base_c(n+10));