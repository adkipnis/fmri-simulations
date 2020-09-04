function w = find_noise_scalar(S_power, N_power, Opts)
    if Opts.sim_type == "mixed"
        % S/N = S/N_emp * w <=> w = (S/N)/(S/N_emp)
        w = Opts.snr ./ (S_power./N_power);
    elseif Opts.sim_type == "signal"
        w = 0;
    elseif Opts.sim_type == "noise"
        w = 1;
    end        
end