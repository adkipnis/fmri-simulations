function p_mat = get_permutation_vectors(Opts)
    perm_idx = zeros(Opts.n_permutations, Opts.n_timepoints-1);
    for i = 1 : Opts.n_permutations
        perm_idx(i,:) = randperm(Opts.n_timepoints-1)+1;
    end
    p_mat = horzcat(ones(Opts.n_permutations, 1), perm_idx);
end