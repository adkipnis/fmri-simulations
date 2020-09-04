function D = generate_sim_4d(S_mat, S_mat_masked, N_mat_masked_scaled, mask_vec, Opts)
    if Opts.sim_type == "noise"
        w = 0;
    else 
        w = 1;
    end
    
    D_mat = S_mat;
    D_mat_masked = w*S_mat_masked + N_mat_masked_scaled;
    D_mat(mask_vec,:) = D_mat_masked;
    D = reshape(D_mat, Opts.size_S);
end