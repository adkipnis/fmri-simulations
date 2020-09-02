function [E, residual_matrix, mask_vector, Opts] = load_residual_matrix(Dirs, Opts)
    [E, nii_hdr] = get_4d_residuals(Dirs);
    Opts.size_E = size(E);
    residual_matrix = reshape(E, [], Opts.size_E(4));
    mask = get_mask(Dirs);
    mask_vector = logical(mask(:));
    Opts.n_voxels = sum(mask_vector);
    Opts.n_timepoints = size(residual_matrix, 2);
    Opts.nii_header = nii_hdr;
end