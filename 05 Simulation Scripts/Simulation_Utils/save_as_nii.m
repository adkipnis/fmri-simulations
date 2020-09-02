function save_as_nii(A, Dirs, Opts, i, r, type)
    for p = 1 : size(A, 5)
        nii = struct('hdr', Opts.nii_header, 'img', A(:,:,:,:,p));
        nii.hdr.hist.descrip = ['sub-', sprintf('%02s', string(i)), '_run-', ...
            sprintf('%02s', string(r)), '_', type ,'_', sprintf('%04s', string(p))];
        nii.hdr.dime.dim(5) = size(A, 4);
        filename = fullfile(Dirs.output_dir, [type, '_', sprintf('%04s', string(p)), '.nii']);
        save_nii(nii, filename)
        fprintf('Saved array to "%s"\n', filename)
    end 
end