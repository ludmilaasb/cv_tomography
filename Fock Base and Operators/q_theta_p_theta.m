% Observáveis q_\theta e p_\theta


theta = 0:1:180;

q_thetas = cell(length(theta),1);

% ---------- dimension fock ---------------
df = 20;
a = annihilfock_mobral(df);
base_c = basecanon(df);

rhot = base_c{1}*base_c{1}';
% ------------------------------------------

for i=1:length(theta)
            q_thetas{i} = (sqrt(2)/2)*((cosd(theta(i))-i*sind(theta(i)))*a+(cosd(theta(i))+i*sind(theta(i)))*a');
end

            
q_0 = q_thetas{1};