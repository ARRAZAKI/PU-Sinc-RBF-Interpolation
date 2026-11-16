function [Pf,diagPU] = PU_SINC_RBF(d,x,bar_x,m_d,phi,w,f,epsilon,mode,cond_limit,eps_bounds)
% mode = 'loocv' (original fminbnd), 'loocv_grid' (grid+raffinement),
%        'loocv_cond' (LOOCV avec seuil de conditionnement),
%        'loocv_ga' (optimisation GA, MATLAB Global Optimization Toolbox).
% epsilon peut etre [], un scalaire, ou un vecteur par patch (tests de perturbation).
if nargin < 9 || isempty(mode), mode = 'loocv'; end
if nargin < 10 || isempty(cond_limit), cond_limit = 1e12; end
if nargin < 11 || isempty(eps_bounds), eps_bounds = [5,15]; end


tilde_x = MakeSDGrid(d,m_d); % Create m_d^d equally spaced PU centres
delta = sqrt(2)/m_d;  supp = 1/delta; % PU radius and parameter for weights
m = size(tilde_x,1); 
s = size(bar_x,1); 
Pf = zeros(s,1);  % Initialize global fit

eps_vec = nan(m,1);
cond_vec = nan(m,1);
cond_eff_vec = nan(m,1);
fallback_vec = false(m,1);

% Shepard matrix
DM_eval = DistanceMatrix(bar_x,tilde_x); 
W = w(supp,DM_eval);  
den = W*ones(m,1);
den(den == 0) = 1; 
W = spdiags(1./den,0,s,s)*W; 

q = ceil(1./delta); % integer-based partitioning parameter

% Partition structures
X_block     = IntegerBasedStructure(x,q,delta,d); 
bar_X_block = IntegerBasedStructure(bar_x,q,delta,d);

% ------------------------------------------------------------
% OPTION: epsilon par patch (si epsilon = [])
% ------------------------------------------------------------
if isempty(epsilon)
    for j = 1:m
        k = IntegerBasedContainingQuery(tilde_x(j,:),q,delta,d);
        [X_N, idx_X_N] = IntegerBasedNeighbourhood(x,X_block,k,q,d);
        n_j = IntegerBasedRangeSearch(tilde_x(j,:),delta,X_N,idx_X_N);
        if isempty(n_j) || numel(n_j) < 3
            eps_vec(j) = mean(eps_bounds);
            fallback_vec(j) = true;
            continue
        end
        [eps_vec(j),fallback_vec(j)] = select_eps_patch_revision( ...
            x(n_j,:),f(n_j),phi,mode,cond_limit,eps_bounds);
    end
elseif isscalar(epsilon)
    eps_vec = epsilon*ones(m,1);
elseif numel(epsilon)==m
    eps_vec = epsilon(:); % garder les memes epsilon pour la perturbation
else
    error('epsilon doit etre [], un scalaire ou un vecteur de m valeurs.');
end

% Ne pas afficher une ligne par patch dans les experiences repetees.

% ------------------------------------------------------------
% Boucle PU (modif minimale : epsilon -> eps_vec(j))
% ------------------------------------------------------------
for j = 1:m % Loop over subdomains

    k = IntegerBasedContainingQuery(tilde_x(j,:),q,delta,d);

    % interpolation points in subdomain
    [X_NeigBlock, idx_X_NeigBlock] = IntegerBasedNeighbourhood(x,X_block,k,q,d);
    n_j = IntegerBasedRangeSearch(tilde_x(j,:),delta,X_NeigBlock,idx_X_NeigBlock);

    % evaluation points in subdomain
    [bar_X_NeigBlock, idx_bar_X_NeigBlock] = IntegerBasedNeighbourhood(bar_x,bar_X_block,k,q,d);
    s_j = IntegerBasedRangeSearch(tilde_x(j,:),delta,bar_X_NeigBlock,idx_bar_X_NeigBlock);

    if (~isempty(s_j)) && (~isempty(n_j))

        epsj = eps_vec(j);

        K_j = phi(epsj, x(n_j,:), x(n_j,:));
        cond_vec(j) = cond(K_j);  
        %lam = 1e-12*max(1,trace(K_j)/size(K_j,1));
        %K_eff = K_j + lam*eye(size(K_j,1));% regularisation
        %cond_eff_vec(j) = 1/rcond(K_eff);
        c_j = K_j\f(n_j); % 

        K_eval = phi(epsj, bar_x(s_j,:), x(n_j,:));
        P_fj = K_eval*c_j;

        Pf(s_j) = Pf(s_j) + P_fj.*W(s_j,j);
    end
end

diagPU.eps_patch  = eps_vec;
diagPU.cond_patch = cond_vec;           % cond(K_j), avant regularisation
diagPU.cond_effective = cond_eff_vec;   % cond(K_j+lambda I)
diagPU.fallback_patch = fallback_vec;    % aucun candidat sous le seuil
idx_used = isfinite(cond_vec);
diagPU.unstable_fraction = mean(cond_vec(idx_used)>cond_limit);
diagPU.centers = tilde_x;

end


function [ep,used_fallback] = select_eps_patch_revision(X,rhs,phi,mode,cond_limit,bounds)
used_fallback = false;
obj = @(v) CostEps_revision(v,X,phi,rhs,mode,cond_limit);
if strcmpi(mode,'loocv_ga')
    if exist('ga','file')~=2
        error('Mode loocv_ga : Global Optimization Toolbox requis (fonction ga).');
    end
    go=optimoptions('ga','PopulationSize',16,'MaxGenerations',12,'Display','off');
    ep=ga(@(v)obj(v(1)),1,[],[],[],[],bounds(1),bounds(2),[],go);
    return
end
if strcmpi(mode,'loocv')
    ep = fminbnd(obj,bounds(1),bounds(2));
    return
end
% Optimisation deterministe par grille + fminbnd dans la meilleure cellule.
gridv = linspace(bounds(1),bounds(2),17);
sc = arrayfun(obj,gridv);
sc(~isfinite(sc)) = Inf;
if all(isinf(sc))
    used_fallback = true;
    obj = @(v) CostEps_revision(v,X,phi,rhs,'loocv',cond_limit);
    sc = arrayfun(obj,gridv); sc(~isfinite(sc))=Inf;
end
[best,idx] = min(sc);
if ~isfinite(best)
    ep = mean(bounds);
    used_fallback = true;
    return
end
lo = gridv(max(1,idx-1)); hi = gridv(min(numel(gridv),idx+1));
if hi>lo
    [eref,score_ref] = fminbnd(obj,lo,hi);
    if isfinite(score_ref) && score_ref<best, ep=eref; else, ep=gridv(idx); end
else
    ep=gridv(idx);
end
end

function J = CostEps_revision(epsilon,X,rbf,rhs,mode,cond_limit)
Kraw = rbf(epsilon,X,X);
n = size(Kraw,1);
if any(~isfinite(Kraw(:)))
    J=Inf; return
end
if strcmpi(mode,'loocv_cond') && rcond(Kraw)<1/cond_limit
    J=Inf; return
end
lam = 1e-12 * max(1,trace(Kraw)/n);
K = Kraw + lam*eye(n);
coef = K\rhs;
[R,p] = chol(K,'lower');
if p==0
    invR = R\eye(n);
    diagInv = sum(invR.^2,1)';
else
    diagInv = diag(K\eye(n));
end
if any(~isfinite(coef)) || any(~isfinite(diagInv)) || any(abs(diagInv)<1e-14)
    J=Inf; return
end
e = coef./diagInv;
J = sqrt(mean(e.^2));
if ~isfinite(J), J=Inf; end
end
