function [N_vec_masked, N] = get_masked_residuals(Dirs, p, mask_vec)
    fname = fullfile(Dirs.input_dir, ['Res_perm_', sprintf('%04s', string(p)), '.nii']);
    N = niftiread(fname);
    N_vec = reshape(N, [], size(N, 4));
    N_vec_masked = N_vec(mask_vec,:);
end