function base_c = basecanon(N)

base_c = cell(N+1, 1); % gera base canonica

for i=1:N+1
    
    base_c {i} = zeros (N+1, 1);
    base_c {i}(i,1) = 1;

end
