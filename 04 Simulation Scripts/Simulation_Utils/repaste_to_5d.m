function D = repaste_to_5d(residual_matrix, y, mask_vector, Opts)    
    new_residual_matrix = repmat(residual_matrix, [1 1 Opts.n_permutations]);
    if size(y,3) ~= size(new_residual_matrix,3)
        y = repmat(y, [1 1 Opts.n_permutations]);
    end
    y = permute(y, [2 1 3]);
%     dif = new_residual_matrix(mask_vector,:,:) - y_perm;
    new_residual_matrix(mask_vector,:,:) = y;
    D = reshape(new_residual_matrix, [Opts.size_E, Opts.n_permutations]);
end
