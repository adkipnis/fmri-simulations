function [y_perm, eps_perm] = permute_residuals(y_pred, eps, p_mat, Opts)
    eps_perm = zeros([size(eps) Opts.n_permutations]);
    for i = 1 : Opts.n_permutations
        eps_perm(:,:,i) = eps(p_mat(i,:)',:);
    end
    y_perm = repmat(y_pred, [1 1 Opts.n_permutations]) + eps_perm;
end