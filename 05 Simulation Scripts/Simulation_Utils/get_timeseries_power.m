function power = get_timeseries_power(A, Opts)
if length(size(A)) > 2
    A_vec = reshape(A, [], size(A, 4));
else
    A_vec = A;
end


power = nanmean(A_vec.^2, 2);
if Opts.snr_type == 'total'
   power = nanmean(power);   
end

end