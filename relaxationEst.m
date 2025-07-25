function [a, r, g, f] = relaxationEst(y, TE, Nrow, Ncol, lambdaA, lambdaR)
% y: 12 x 4096 
% TE: 12x1
% Nrow: 64
% Ncol: 64

AL_iter = 5e3;
mu = 1e-1;

Npix = size(y, 2); % Npix: 4096
t = 1;
res = inf;
tol = 1e-3;
tol1 = sqrt(Npix)*tol;
tol2 = sqrt(Npix)*tol;
res_p = inf;
res_d = inf;

a = max(y(:)) * rand(Npix, 1); % 4096 x 1
r = 0.1*ones(Npix, 1); % 4096x1

g = a;
f = r;

d1 = zeros(size(g));  % Lagrange Variables   % 4096 x 1
d2 = zeros(size(f));    % 4096 x 1

puxg = 0.*reshape(g, Nrow, Ncol);
puyg = 0.*reshape(g, Nrow, Ncol);

puxr = 0.*reshape(g, Nrow, Ncol);
puyr = 0.*reshape(g, Nrow, Ncol);

disp(' ------- Running T2* Relaxation Estimation ------- ')

while (t<AL_iter) && ((abs (res_p) > tol1) || (abs (res_d) > tol2))
    g01 = g;
    f02 = f;
    
    %% Solve for g
    % %%%%%% INSERT YOUR CODE HERE % %%%%%%
    varargin_g = lambdaA/mu;
    [g_img,dual_x,dual_y,err] = chambolle_prox_TV_stop(reshape(a - d1, Nrow, Ncol), varargin_g);
    g = reshape(g_img, [], 1);

    
    %% Solve for f
    % %%%%%% INSERT YOUR CODE HERE % %%%%%%
    varargin_f = lambdaR/mu;
    [f_img,~,~,~] = chambolle_prox_TV_stop(reshape(r - d2, Nrow, Ncol), varargin_f);
    f = reshape(f_img, [], 1);
    
    for oo = 1:Npix
            %% Solve for a
            % %%%%%% INSERT YOUR CODE HERE % %%%%%%
            numer = ((y(:, oo)' * exp(-r(oo) * TE)) + mu * (g(oo) + d1(oo)));
            denom = (exp(-r(oo) * TE)' * exp(-r(oo) * TE) + mu);
            a(oo) = numer / denom;
            

        %% Solve for r
            % %%%%%% INSERT YOUR CODE IN THE gradientDescentAK FUNCTION % %%%%%%
            %[r(oo), niter, rnorm, dx] = gradientDescentAK(r(oo), TE, a(oo), y(:, oo), mu, f_d2(oo));

            f_d2 = f(oo) - d2(oo);
            [r(oo), niter, rnorm, dx] = gradientDescentAK(r(oo), TE, a(oo), y(:, oo), mu, f_d2);
    end
   
    %% Update Lagrange multipliers
    %disp('lagrange update step')
    %disp(size(a)+' '+size(g));
    res1 = a - g;
    d1 = d1 - res1;
    
    res2 = r - f;
    d2 = d2 - res2;
    
    %% Update mu
    if mod(t, 10) == 1
        res_p = sqrt(norm(res1, 'fro')^2 + norm(res2, 'fro')^2);    % primal residue
        res_d = mu*(norm(g01 - g + f02 - f,'fro'));                 % dual residue
        
        % update mu
        if res_p > 10*res_d
            mu = mu*2;
            d1 = d1/2;
            d2 = d2/2;
        elseif res_d > 10*res_p
            mu = mu/2;
            d1 = d1*2;
            d2 = d2*2;
        end
    end
    res(1) = norm(res1, 'fro');
    res(2) = norm(res2, 'fro');
    
    if mod(t, 10) == 1
        fprintf(' t = %f, res_p = %f, res_d = %f, mu = %f, niterGD = %f \n', t, res_p, res_d, mu, niter)
    end
    
    t = t + 1;
    
end
disp(' ------- T2* Relaxation Estimation Completed ------- ')

end