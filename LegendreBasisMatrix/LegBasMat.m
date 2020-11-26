clearvars;
%%
T_Leg = tic;
N = 1024;
N_vec0 = (0 : N);
N_basis = N / 8;
fprintf('Generating symbolic Legendre functions for n = 0 : %d and basis functions for n = 0 : %d.\n', N * 2, N_basis);    

x = sym('x');
P  = cell([1,2*N+1]);
Q  = cell([1,2*N_basis+1]);
W  = cell([1,2*N_basis+1]);
Pi = cell([1,2*N_basis+1]);
Qi = cell([1,2*N_basis+1]);
Wi = cell([1,2*N_basis+1]);

T_n = tic;
n = 0;
fprintf('\tGenerating functions for n = %04d. ', n); 
P{ indP(n)} = sym(1);
W{ indP(n)} = sym(0);
Q{ indP(n)} = 1/2*log((1+x)/(1-x));
Pi{indP(n)} = sym(1);
Wi{indP(n)} = sym(0);
Qi{indP(n)} = -1i * atan(x);
fprintf('Time: %s.\n',datestr(seconds(toc(T_n)),'SS.FFF'));

T_n = tic;
n = 1;  
fprintf('\tGenerating functions for n = %04d. ', n);
P{ indP(n)} = x;
W{ indP(n)} = sym(1);
Q{ indP(n)} = P{indP(n)} * Q{indP(0)} - W{indP(n)};
Pi{indP(n)} = x/1i;
Wi{indP(n)} = sym(1);
Qi{indP(n)} = Pi{indP(n)} * Qi{indP(0)} - Wi{indP(n)};
fprintf('Time: %s.\n',datestr(seconds(toc(T_n)),'SS.FFF'));

for n = 1 : N*2-1
    T_n = tic;
    fprintf('\tGenerating functions for n = %04d. ', n+1); 
        P{ indP(n+1)} = expand( (2*n+1)*  x     * P{ indP(n)} - n * P{ indP(n-1)}) / (n+1);
    if n <= N_basis*2-1
        W{ indP(n+1)} = expand( (2*n+1)*  x     * W{ indP(n)} - n * W{ indP(n-1)}) / (n+1);
        Q{ indP(n+1)} = P{ indP(n+1)} * Q{ indP(0)} - W{ indP(n+1)};
        
        Pi{indP(n+1)} = expand( (2*n+1)* (x/1i) * Pi{indP(n)} - n * Pi{indP(n-1)}) / (n+1);
        Wi{indP(n+1)} = expand( (2*n+1)* (x/1i) * Wi{indP(n)} - n * Wi{indP(n-1)}) / (n+1);
        Qi{indP(n+1)} = Pi{indP(n+1)} * Qi{indP(0)} - Wi{indP(n+1)};
    end
    fprintf('Time: %s.\n',datestr(seconds(toc(T_n)),'SS.FFF'));
end

PP  = P( indP(0:2:N*2));
P   = P( indP(0:N_basis*2));
Q   = Q( indP(0:N_basis*2));
PPi = Pi(indP(0:2:N_basis*2));
QQi = Qi(indP(0:2:N_basis*2));
fprintf('Time for generating Legendre functions: %s.\n',datestr(seconds(toc(T_Leg)),'HH:MM:SS.FFF'));
%%
fprintf('Saving Legendre functions. ');
save('LegendreFunctions.mat','N','N_basis','x','P','Q','PP');
fprintf('Total time for Legendre functions: %s.\n\n',datestr(seconds(toc(T_Leg)),'HH:MM:SS.FFF'));

%%
T_A = tic;
if ~exist(sprintf('A_partial_%d',N),'dir')
    mkdir(sprintf('A_partial_%d',N));
end
if ~exist('out_err','dir')
    mkdir('out_err');
end
if ~exist('PCstorage','dir')
    mkdir('PCstorage');
end
fprintf('Submitting jobs to calculate A matrix.\n');
filename = 'A0_matrix.slurm';
write_slurm_file(filename,N);
system(sprintf('sbatch %s',filename));
fprintf('\n');
%%
T_B = tic;
fprintf('Calculating symbolic coefficients and basis functions for n = 0 : %d.\n',  N_basis);

syms eta xi

CP = cell([1,N_basis+1]);
CQ = cell([1,N_basis+1]);
Mdiff = cell([1,N_basis+1]);
MM = cell([1,N_basis+1]);

U0  = cell([1,N_basis+1]);
J0  = cell([1,N_basis+1]);
U   = cell([1,N_basis+1]);

for n = 0 : N_basis
    T_N = tic;
    fprintf('\tCalculating n = %04d. ',n);
    U0{indP(n)}= subs(PP{indP(n)}, x, eta);
    
    if n ==0 
        CP{indP(n)} = sym(1);
    else
        CP{indP(n)} = CP{indP(n-1)} * (-1) / (2 * n - 1) * (2 * n);
    end
    CQ{indP(n)} = CP{indP(n)} * (sym(-2i)/sym(pi));
    MM{indP(n)}=  simplify(CP{indP(n)} * PPi{indP(n)} + CQ{indP(n)} * QQi{indP(n)});
    U{indP(n)} = (subs(MM{indP(n)}, x, xi) * subs(PP{indP(n)}, x, eta));
    
    Mdiff{indP(n)} = sym(-1i) * CP{indP(n)} * CQ{indP(n)};
    J0{indP(n)}= (-Mdiff{indP(n)} * subs( PP{indP(n)} / x, x, eta)); 
    fprintf('Time: %s.\n',datestr(seconds(toc(T_N)),'SS.FFF'));
end
fprintf('Time for symbolic coefficients and basis functions: %s.\n', datestr(seconds(toc(T_B)),'MM:SS.FFF'));
%%
fprintf('Saving basis functions. ');
save('BasisFunctions.mat', 'N_basis', 'N_vec0', 'eta', 'xi', 'CP', 'U0', 'J0', 'U');
fprintf('Total time for basis functions: %s.\n\n',datestr(seconds(toc(T_B)),'MM:SS.FFF'));

%%
T_Mat = tic;
fprintf('Calculating numeric CP, CQ, M''(0) and M matrices.\n');

CP_sym = (cellfun(@eval,CP));
CQ_sym = -2i/pi * CP_sym;
Mdiff_sym = -1i * CP_sym .* CQ_sym;

CP_num = zeros(1,N+1);
for n = N_vec0
    if n ==0 
        CP_num(indP(n)) = (1);
    else
        CP_num(indP(n)) = CP_num(indP(n-1)) * (-1) / (2 * n - 1) * (2 * n);
    end
end
CQ_num = CP_num * (-2i/pi);
Mdiff_num = -1i * CP_num .* CQ_num;

M0 = diag( - Mdiff_num * pi/4 ./ (4 * N_vec0 + 1 ) );
% M1 = M0(indP(1):indP(N),indP(1):indP(N));

fprintf('Saving CP, CQ, M''(0) and M matrices. ');
save('MatrixData.mat','N','N_basis','N_vec0','CP_num','CQ_num','Mdiff_num','CP_sym','CQ_sym','Mdiff_sym','M0');
fprintf('Total time : %s.\n\n',datestr(seconds(toc(T_Mat)),'SS.FFF'));

%%
fprintf('Check status of A sub-matrices calculation.\n');
D = dir(fullfile(sprintf('A_partial_%d',N),'*.mat'));

while length(D) < (N + 1)
    fprintf('%d .mat files found. %d files sub-matrices calculation ongoing...', length(D), (N+1)-length(D));
    fprintf('Time: %s. Check status again in 5 minutes...\n', datestr(seconds(toc(T_A)),'HH:MM:SS.FFF'));
    pause(300);
    D = dir(fullfile(sprintf('A_partial_%d',N),'*.mat'));
end

fprintf('%d .mat files found. A sub-matrices calculation finished. Time: %s.\n\n', length(D), datestr(seconds(toc(T_A)),'HH:MM:SS.FFF'));

%%
T_Mat = tic;
fprintf('Assembling A matrix from sub-matrices.\n');
A0 = zeros(N+1,N+1);

for n = 0 : N
    T_Ap = tic;
    fprintf('\tCalculating n = %04d.',n);
    load(fullfile(sprintf('A_partial_%d',N),sprintf('A_%d.mat', n)),'A0_p');
    A0 = A0 + A0_p;
    fprintf(' Time: %s.\n',datestr(seconds(toc(T_Ap)),'MM:SS.FFF'));
end  
fprintf('Time for assembling matrix A: %s.\n', datestr(seconds(toc(T_Mat)),'HH:MM:SS.FFF'));

% A1 =  A0(indP(1):indP(N),indP(1):indP(N));
% a0  = A0(indP(0):indP(N),indP(0));
% a1  = A0(indP(1):indP(N),indP(0));

fprintf('Saving A matrix.\n');
save('MatrixData.mat','A0','-append');
fprintf('Total time for matrix A: %s.\n',datestr(seconds(toc(T_A)),'HH:MM:SS.FFF'));
