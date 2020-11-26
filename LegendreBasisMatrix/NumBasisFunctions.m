if ~exist('N_basis ','var') || ~exist('eta','var') || ~exist('xi','var') || ~exist('U','var')|| ~exist('U0','var')|| ~exist('J0','var')
    T_M = tic;
    fprintf('Loading Basis Functions. ');
    load('BasisFunctions.mat', 'N_basis', 'eta', 'xi', 'U0', 'J0', 'U');
    fprintf('Time: %s.\n\n',datestr(seconds(toc(T_M)),'MM:SS.FFF'));  
end

%%
T_B = tic;
N_plot = 10;
digitsOld = digits(128);
fprintf('Calculating numeric basis functions U0 and J0 for n = 0 : %d.\n',  N_basis );

delta_r = 0.0001;
r_num = [0   : delta_r * 10 : 0.99 - delta_r * 10, ...
         0.99 : delta_r  : 1 ]';
z_num = 0;
U0_num = zeros(length(r_num), N_basis +1);
J0_num = zeros(length(r_num)-1, N_basis +1);

eta_num = sqrt( ( sqrt( ((r_num - 1).^2 + z_num^2).*((r_num + 1).^2 + z_num^2 ) ) - (r_num.^2 + z_num^2 - 1) ) / 2);

for ii = 0 : N_basis  
    T_N = tic;
    fprintf('\tCalculating n = %02d. ',ii);
    U0_num(:,indP(ii)) = subs( U0{indP(ii)}, eta, eta_num);
    J0_num(:,indP(ii)) = subs( J0{indP(ii)}, eta, eta_num(1:end-1));
    fprintf('Time: %s.\n',datestr(seconds(toc(T_N)),'MM:SS.FFF'));
end
fprintf('Time for numeric basis functions U0 and J0: %s.\n', datestr(seconds(toc(T_B)),'HH:MM:SS.FFF'));
fprintf('Saving basis functions U0 and J0. ');
save('BasisFunctions_U0J0_num.mat', 'N_basis' ,'N_plot', 'r_num', 'J0_num', 'U0_num');
fprintf('Total time for basis functions U0 and J0: %s.\n\n',datestr(seconds(toc(T_B)),'HH:MM:SS.FFF'));


%%
T_B = tic;
N_plot = (0:9);
fprintf('Calculating numeric basis functions U for n = 0 : %d.\n',  N_plot(end));

r_num_end = [8, 3.5, 3.5, 2, 2, 2, 1.25, 1.25, 1.25, 1.25];
z_num_end = r_num_end;


R_num = cell(size(N_plot));
Z_num = cell(size(N_plot));
U_num = cell(size(N_plot));

delta_r = 0.0001;
delta_z = 0.0001;
eps_num = 1e-14;

r_num_fine = [  0                   : delta_r * 100 : 0.9, ...          % 0.01, 91
                0.9 + delta_r * 10  : delta_r * 10 : 0.99, ...          % 0.001, 90
                0.99 + delta_r * 1  : delta_r * 1 :  1, ...             % 0.0001, 100
                1.0 + delta_r * 10  : delta_r * 10 :  1.01 ];          % 0.001, 10
z_num_fine = [  0     : delta_z * 1       : 0.01 - delta_z * 1, ...    % 0.0001, 100
                0.01  : delta_z * 10      : 0.1 - delta_z * 10, ...    % 0.001, 90
                0.1   : delta_z * 100     : 1  - delta_z * 100, ...    % 0.01, 90
                1   ];                 % 
%

for ii = 1 : length(N_plot) 
    T_N = tic;
    fprintf('\tCalculating n = %02d.\n',N_plot(ii));
    
    if (ii == 1) || ( (ii > 1) &&  (r_num_end(ii) ~= r_num_end (ii-1) )  )
        fprintf('\t\tGenerating new sampling coordinates. ');
        r_num = [   r_num_fine, ...
                    1.01 + delta_r * 100 : delta_r * 100 : min(2, r_num_end(ii)), ...
                    min(2, r_num_end(ii)) + delta_r * 1000 : delta_r * 1000 : r_num_end(ii) ]' ;    % 0.1, 79    
        z_num = [   z_num_fine, ...
                    1 + delta_z * 1000 : delta_z * 1000 : z_num_end(ii) ];


        [R_num{ii}, Z_num{ii}] = ndgrid(r_num, z_num);
        tmp_1 = sqrt( ((1 - R_num{ii}).^2 + Z_num{ii}.^2) .* ((R_num{ii} + 1).^2 + Z_num{ii}.^2 ) );
        tmp_2 = (R_num{ii}.^2 + Z_num{ii}.^2 - 1);
        tmp_3 = ( tmp_1 + tmp_2 ) / 2; tmp_3(abs(tmp_3) < eps_num) = 0;
        tmp_4 = ( tmp_1 - tmp_2 ) / 2; tmp_4(abs(tmp_4) < eps_num) = 0;

        Xi_num  = (sqrt( tmp_3 ));
        Eta_num = (sqrt( tmp_4 ));
        fprintf('r: %d, z: %d. \n', length(r_num), length(z_num));
    else
        R_num{ii} = R_num{ii-1};
        Z_num{ii} = Z_num{ii-1};
    end
    fprintf('\t\tEvaluating basis function.');
    fac = factor(U{indP(N_plot(ii))},[eta, xi]);
    
    if ii ~= 1
        U_num{ii} = double(vpa(fac(1))) .* double(vpa( subs( fac(2), eta, Eta_num) )) .* ...
                                           double(vpa( subs( fac(3), xi,  Xi_num ) )) ;
    else
        U_num{ii} = double(vpa(fac(1))) .* double(vpa( subs( fac(2), xi, Xi_num) ));
    end
    fprintf('Time: %s.\n',datestr(seconds(toc(T_N)),'MM:SS.FFF'));
%     U_num{ii}(:,:) = eval( subs( U{indP(N_plot(ii))}, {xi, eta}, {Xi_num, Eta_num}) );
end
digits(digitsOld)
fprintf('Time for numeric basis functions U: %s.\n', datestr(seconds(toc(T_B)),'MM:SS.FFF'));

fprintf('Saving basis functions U. ');
save('BasisFunctions_U_num.mat','N_plot','R_num','Z_num','U_num','r_num_end','z_num_end');
fprintf('Total time for basis functions U: %s.\n\n',datestr(seconds(toc(T_B)),'MM:SS.FFF'));
