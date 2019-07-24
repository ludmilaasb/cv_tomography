function [base_fock_bosons] = base_Fock_bosons_multimod(N,dbase)


%Incompleto! Não usar!
dim_fock = dbase^N - nchoosek(dbase,N);

base_c = cell(dim_fock, 1); % gera base canonica

for i=1:dim_fock
    
    base_c {i} = zeros (dim_fock, 1);
    base_c {i}(i,1) = 1;

end

modctrl = zeros(N,1);
for i=1:N
    modctrl(i) = dbase;
end


base_fock_bosons = cell(modctrl);
k=1;

% for aloop = 1:numel(base_fock_bosons)
%    for ii = 1:dbase
%        if ii>=ii
%                 base_fock_bosons{i,j} = base_c{k};                
%                 k = k+1;
%                 base_fock_bosons{j,i} = base_fock_bosons{i,j};
%        end
%    end
% end
% 
% 
%    
%     for i=1:dbase
%         for j=1:dbase
%             if j>=i
%                 base_fock_bosons{i,j} = base_c{k};                
%                 k = k+1;
%                 base_fock_bosons{j,i} = base_fock_bosons{i,j};
%             end;
%         end;
%     end;
% end;