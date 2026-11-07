function ind_out = unique_fix(ind0, N)


    ind = round(ind0(:));
    M = numel(ind);
    ind_out = ind;
    
    ind_out(ind_out < 1) = 1;
    ind_out(ind_out > N) = N;

    [~, first_occurrence] = unique(ind_out, 'stable');
    duplicate_idx = setdiff(1:M, first_occurrence);

    used = ind_out(first_occurrence);
    available = setdiff(1:N, used, 'stable');
    
    for i = 1:length(duplicate_idx)
        if isempty(available)
            break;
        end
        ind_out(duplicate_idx(i)) = available(1);
        available(1) = [];
    end
    
    if isrow(ind0)
        ind_out = ind_out.';
    end
end
