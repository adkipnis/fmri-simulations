function p_mat = get_permutation_vectors(Opts)
    perm_idx = zeros(Opts.n_permutations, Opts.n_timepoints-Opts.ar_n);
    for i = 1 : Opts.n_permutations
        perm_idx(i,:) = randperm(Opts.n_timepoints-Opts.ar_n)+Opts.ar_n;
    end
    p_mat = horzcat(repmat([1:Opts.ar_n], Opts.n_permutations, 1), perm_idx);
end