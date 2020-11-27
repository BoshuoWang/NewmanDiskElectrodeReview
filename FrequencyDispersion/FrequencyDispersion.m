%%
close all;
addpath(fullfile('..','LegendreBasisMatrix'));

if ~exist('N','var') || ~exist('A0','var') || ~exist('M0','var')
    T_A = tic;
    fprintf('Loading matrices.');
    load('MatrixData.mat','N','A0','M0');
    fprintf('Time: %s.\n\n',datestr(seconds(toc(T_A)),'MM:SS.FFF'));
end

%%
T_F = tic;

log_range = 5;
num_per_dec = 12;

No_N = round(2.^(-6:0) * N);

Omega = [0, logspace(-log_range,log_range, log_range*2*num_per_dec+1 )];
G = [0, 0.01, 0.1, 1, 10];

N_omega = length(Omega);
N_G = length(G);

bH = cell(size(No_N));

BH0  = zeros(N_omega, N_G, length(No_N));
Reff = zeros(N_omega, N_G, length(No_N));
Ceff = zeros(N_omega, N_G, length(No_N));

a0  = A0(indP(0):indP(N),indP(0));

for kk =  length(No_N): -1 : 1
    T_N = tic;
    fprintf('Calculating frequency dispersion with N = %d.\n', No_N(kk));
    ind = indP( 0 : No_N(kk) );
    A0_N = A0( ind , ind );
    M0_N = M0( ind , ind );
    a0_N = a0(ind);
    
    bH{kk} = zeros( indP(No_N(kk)), N_omega, N_G);
    for ii = 1 : N_G
        fprintf('\tCalculating G = %2.2f.\n', G(ii));
        
        for jj = 1 : N_omega
            bH{kk}(:, jj, ii) = ( A0_N * (G(ii)+1i*Omega(jj)) + M0_N ) \ ((G(ii) + 1i * Omega(jj)) * a0_N) ;
            
            BH0(jj,ii,kk) = bH{kk}(indP(0),jj,ii);
            Reff(jj, ii, kk) = real(BH0(jj,ii,kk))./abs(BH0(jj,ii,kk)).^2;
            Ceff(jj, ii, kk) = imag(BH0(jj,ii,kk)) /abs(BH0(jj,ii,kk)).^2 * Omega(jj);
        end
        
    end
    
    Reff(1,1,kk) = 32/(3*pi^2);
    Ceff(1,1,kk) = 1;

    if kk ~= length(No_N)
        fprintf('\tRelative errors of BH0 vs. N = %d: \n',N);
        err = squeeze(abs( BH0(:,:,kk)./BH0(:,:,end) -1 ));
        fprintf('\tBH0: Max: %2.2g. Min: %2.2g. Mean: %2.2g. StdDev: %2.2g.\n', ...
            max(err(:)), min(err(:)), nanmean(err(:)),  nanstd(err(:)));
        err = squeeze(abs( bH{kk}(:,:,:)./bH{end}(ind,:,:) -1 ));
        fprintf('\tbH:  Max: %2.2g. Min: %2.2g. Mean: %2.2g. StdDev: %2.2g.\n', ...
            max(err(:)), min(err(:)), nanmean(err(:)),  nanstd(err(:)));
    end
 
    fprintf('Time: %s.\n\n',datestr(seconds(toc(T_N)),'SS.FFF'));
    
end


BH_inf = A0 \ a0;
Reff_inf = real(BH_inf(1))./abs(BH_inf(1)).^2;
Ceff_inf = inf;

fprintf('Saving frequency dispesion results. ');
save('freq_disp_data.mat','G','Omega','bH','BH0','Reff','Ceff','Reff_inf','Ceff_inf');
fprintf('Total time for frequency dispesion: %s.\n\n',datestr(seconds(toc(T_F)),'MM:SS.FFF'));

%%
rmpath(fullfile('..','LegendreBasisMatrix'))
