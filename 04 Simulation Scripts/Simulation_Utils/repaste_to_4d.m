function D = repaste_to_4d(residual_matrix, y, mask_vector, Opts)    
    residual_matrix(mask_vector,:,:) = y';
    D = reshape(residual_matrix, Opts.size_E);
end
