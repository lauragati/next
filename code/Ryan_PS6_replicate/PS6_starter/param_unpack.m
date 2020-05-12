param_list = fieldnames(param);
param_val = struct2cell(param);
idx_yes = zeros(1,length(param_list));
for j = 1:length(param_list)
    idx_yes(j) = sum(~isnan(param_val{j}));
end


for j = find(idx_yes==1)
    eval([param_list{j} '= param_val{j};']);
    %param_list{j} = param_val{j};
end

if exist('set', 'var')
    
    set_list = fieldnames(set);
    set_val = struct2cell(set);
    idx_yes = ones(1,length(set_list));
    for j = 1:length(set_list)
     %   idx_yes(j) = sum(~isnan(set_val{j}));
    end
    
    for j = find(idx_yes>0)
        eval([set_list{j} '= set_val{j};']);
    end
    
end