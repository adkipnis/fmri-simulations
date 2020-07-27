function basis_multi = effects_of_interest_multi(Opts, n)
    
    basis = zeros(Opts.n_reg, Opts.n_reg-1);
    step = 1 + sum(Opts.hrf_derivs);

    for i = 1:step:Opts.n_stim_betas
        basis(i,i) = 1;
    end
    
    basis_multi_tmp = repmat(basis, 1, n);
    basis_multi = [basis_multi_tmp, zeros(Opts.n_reg, n)];
    
end