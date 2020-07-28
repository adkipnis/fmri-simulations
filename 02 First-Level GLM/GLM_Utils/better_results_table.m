function better_results_table(TabDat, results_dir, contrast_name,  Opts)
    % Partially derived from https://github.com/solleo/myspm/edit/master/myspm_result.m

    % read SPM.mat to find number of regressors and number of contrast matrices
    alpha = Opts.alpha_level;
    table_path =  [results_dir,'/TabDat_', contrast_name,'.csv']; 
    footer_path =  [results_dir,'/TabDat_', contrast_name,'_footer.txt']; 

    % table headers
    fid = fopen(table_path,'w'); % for other applications
    fprintf(fid, cell2fmt(TabDat.hdr(1,:)));
    fprintf(fid, strrep(strrep(strrep(cell2fmt(TabDat.hdr(3,:)),'\it',''),'\rm_',''),'\equiv',''));

    % gather data
    n_peaks = size(TabDat.dat,1);
    COORDS={};
    PI=[];
    pvals=[];
    tmax=[];
    zmax=[];

    for peak_index = 1:n_peaks
        pvals = [pvals TabDat.dat{peak_index,7}];
        tmax = [tmax TabDat.dat{peak_index,9}];
        zmax = [zmax TabDat.dat{peak_index,10}];
        n_cols = size(TabDat.dat,2);
        for col = 1:n_cols
          fprintf(fid, TabDat.fmt{col},TabDat.dat{peak_index,col});
          if col<n_cols
            fprintf(fid, '\t');
          else
            fprintf(fid, '\n');
          end
        end
        if (TabDat.dat{peak_index,7} < alpha)
          COORDS = [COORDS TabDat.dat{peak_index,end}];
          PI = [PI peak_index];
        end
    end 
    fclose(fid);

    % Footer
    fid = fopen(footer_path,'w'); 
    for footnote = 1:length(TabDat.ftr)
        fprintf(fid, cell2fmt(TabDat.ftr(footnote,1)), TabDat.ftr{footnote,2});
    end
    fclose(fid);


    fprintf('Printed results table to "%s".\n', table_path)
    fprintf('...and its footer to "%s".\n', footer_path)
    
end