function filename = save_as_nii(A, Dirs, Opts, i, r, p, type)
    if length(size(A)) < 5
            nii = struct('hdr', Opts.nii_header, 'img', A);
            nii.hdr.hist.descrip = ['sub-', sprintf('%02s', string(i)), '_run-', ...
                sprintf('%02s', string(r)), '_', type ,'_', sprintf('%04s', string(p))];
            nii.hdr.dime.dim(1) = length(size(A));
            nii.hdr.dime.dim(5) = size(A, 4);
            filename = fullfile(Dirs.output_dir, [type, '_', sprintf('%04s', string(p)), '.nii']);
            save_nii(nii, filename)
            fprintf('Saved array to "%s"\n', filename)
    elseif length(size(A)) == 5
        for p = 1 : size(A, 5)
            A_tmp = A(:,:,:,:,p);
            save_as_nii(A_tmp, Dirs, Opts, i, r, p, type);
        end
    else
        warning('Input data has too many dimensions. No nii will be written.')
    end

end