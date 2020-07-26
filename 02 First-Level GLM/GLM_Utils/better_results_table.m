function better_results_table(Dirs, spm_results, TabDat)
    % Partially derived from https://github.com/solleo/myspm/edit/master/myspm_result.m

%     % set table filename
%     today = datestr(now,'yyyymmmdd');
%     csv_path = fullfile(Dirs.results_dir,['spm_',today,'.csv']);
%     fid = fopen(csv_path,'w');
%     headers = 'measure\tcontrast_name\teffect_size\tstat\tpeak\tpeakZ\tuncor_pval\tcor_pval_peak\tcor_pval_clus\tK_E\tx_mm\ty_mm\tz_mm\tpeak_strc_name\tpeak_strc_prob\n';
%     fprintf(fid, headers, spm_results.matlabbatch{1}.spm.stats.results.conspec.threshdesc);
%     fclose(fid);

    % read SPM.mat to find number of regressors and number of contrast matrices
    alpha = spm_results.matlabbatch{1}.spm.stats.results.conspec.thresh;
    load([Dirs.results_dir,'/SPM.mat']);
%     n_reg = size(SPM.xX.X,2);
    n_con = numel(SPM.xCon);   
    
    for i=1:n_con
        con_title{i} = SPM.xCon(i).name; % Contrast analysis title
        table_path =  [Dirs.results_dir,'/TabDat_',con_title{i},'.csv']; 
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
        fprintf('Printed table to "%s".\n', table_path)
    end
end