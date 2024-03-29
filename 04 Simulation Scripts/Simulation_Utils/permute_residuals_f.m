function y_perm = permute_residuals(A, eps, p_mat, Opts)
% Permute AR residuals (eps) and reuse AR model parameters per voxel to
% generate "permuted" GLM residuals
    y_perm = zeros([size(eps) Opts.n_permutations]);
    for i = 1 : Opts.n_permutations
        eps_perm = eps(p_mat(i,:)',:);
        for k = 1 : size(eps,2)
            y_perm(:,k,i) = filter(1,[1 A(:,:,k)], eps_perm(:,k)); % altered by HS
        end
    end
end