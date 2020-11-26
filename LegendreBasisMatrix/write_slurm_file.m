function write_slurm_file(filename,N)

fid = fopen(filename,'w');

fprintf(fid, '#!/bin/bash\n');
fprintf(fid, '#\n');
fprintf(fid, '#SBATCH --array=0-%d%%60\n',N);
fprintf(fid, '#SBATCH --job-name=A_%d\n',N);
fprintf(fid, '#SBATCH --output=out_err/%%a.out\n');
fprintf(fid, '#SBATCH --error=out_err/%%a.err\n');
fprintf(fid, '#SBATCH --mail-type=END,FAIL\n');
fprintf(fid, '#SBATCH --mem=30G\n');        % change parameters according to local cluster 
fprintf(fid, '#SBATCH -c 10\n');            % change parameters according to local cluster
fprintf(fid, newline);
fprintf(fid, 'uname -n 1>&2 \n');
fprintf(fid, newline);
fprintf(fid, '# Local work directory for parfor loop temporary files\n');
fprintf(fid, 'mkdir -p -v PCstorage/$SLURM_JOB_ID\n');
fprintf(fid, 'echo "Created parfor tmp folder"\n');
fprintf(fid, newline);
fprintf(fid, 'n=$(($SLURM_ARRAY_TASK_ID))\n');
fprintf(fid, newline);
fprintf(fid, '/usr/bin/time -v /admin/apps/rhel7/matlabR2019a/bin/matlab -nodisplay -r "partial_A(${n}); quit;"\n');
fprintf(fid, newline);  % change matlab directory and version according to local cluster
fprintf(fid, '# Delete work directory for parfor loop temporary files\n');
fprintf(fid, 'rm -rf -v PCstorage/$SLURM_JOB_ID\n');
fprintf(fid, 'echo "Removed parfor tmp folder"\n');
fprintf(fid, newline);

fclose(fid);
end
