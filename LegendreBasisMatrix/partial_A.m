function A0_p = partial_A(n)
    
    T_Leg = tic;
    fprintf('Loading Legendre Functions. ');
    load('LegendreFunctions','N','x','PP')
    fprintf('Time: %s.\n\n',datestr(seconds(toc(T_Leg)),'MM:SS.FFF'));
    
    pc = parcluster('local');
    numcores = feature('numcores');
    pc_storage_dir = fullfile('PCstorage',getenv('SLURM_JOB_ID'));
    pc.JobStorageLocation =  pc_storage_dir;
    poolobj = parpool(pc, max( 1, min( numcores - 1, N-n) ) );
    fprintf('\n');
    
    T_A = tic;
    fprintf('Calculating partial A matrix: n = %04d.\n',n);
    
    A0_p = zeros(N+1,N+1);
    A0_p(indP(n),indP(n)) = int(PP{indP(n)}* PP{indP(n)}* x, x, 0, 1);

    fprintf('\tm = %04d. AH(n,m) = % 1.4e. Time: %s.\n', n, A0_p(indP(n),indP(n)), datestr(seconds(toc(T_A)),'MM:SS.FFF'));
    temp = zeros(1,N-n);
    
    parfor m = n + 1 : N
        T_m = tic;
        x = sym('x');
        temp(m-n) = int(PP{indP(n)} * PP{indP(m)}* x, x ,0, 1);
        fprintf('\tm = %04d. AH(n,m) = % 1.4e. Time: %s.\n', m, temp(m-n), datestr(seconds(toc(T_m)),'MM:SS.FFF'));
    end

    A0_p(indP(n),indP(n+1:N)) = temp;
    A0_p(indP(n+1:N),indP(n)) = temp';
    fprintf('Total time: %s.\n\n', datestr(seconds(toc(T_A)),'HH:MM:SS.FFF'));

    save(fullfile(sprintf('A_partial_%d',N),sprintf('A_%d.mat', n)),'A0_p');

    if exist('poolobj','var')
        delete(poolobj);
    end
end
