function [N_vec_masked, N] = get_masked_residuals(p, mask_vec)
    N = niftiread(['Res_perm_', sprintf('%04s', string(p)), '.nii']);
    N_vec = reshape(N, [], size(N, 4));
    N_vec_masked = N_vec(mask_vec,:);
end