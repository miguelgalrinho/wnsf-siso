function [default_model,my_model] = wnsfBJ(z,m,varargin)
%wnsfOE - Estimate a SISO Box-Jenkins model with WNSF
%y(t) = B(q)/F(q) u(t) + C(q)/F(q) e(t),
%where B(q) =     b_1 q^{-1} + ... + b_{mb} q^{-mb} 
%      F(q) = 1 + f_1 q^{-1} + ... + f_{mf} q^{-mf}
%      C(q) = 1 + c_1 q^{-1} + ... + c_{mc} q^{-mc}
%      D(q) = 1 + c_1 q^{-1} + ... + d_{md} q^{-md}
%
% Syntax:  [model_idpoly_format,model_alternative_format] = wnsfBJ(data,polynomial_orders,Options)
%
% Inputs:
%    data - Input and Output signals in iddata format
%    polynomial_orders - Vector with two or four entries [mf mb mc md],
%    where mc and md are optional. If mc and md are not provided, the 
%    'UnstructuredOrders' (optional) - Followed by a vector of
%       non-parametric ARX model orders to be tried out (default: 10 linearly
%       spaced orders between 20 and sample_size/4)
%    'Iterations' (optional) - Followed by an integer with the maximum
%       number of WLS iterations (default: 20).
%    'Tol' (optional) - Followed by float with parameter change tolerance
%       for stopping iterations (default: 1e-4)
%    'InitialConditions' (optional) - Folloed by string 'zero', 'truncate',
%       on how to handle initial conditions (default: 'truncate')
%
% Outputs:
%    model_idpoly_format - Estimated model in idpoly format
%    model_alternative_format - Estimated model as a structure with
%    following entries:
%           G - Plant model
%           CostFunction - PEM cost function value for the estimated model
%              (estimated noise variance)
%           BestUnstructuredOrder - Non-parametric ARX model order used for
%              the final parametric estimate
%           LastIteration - Number of weighted least-squares iterations
%              performed for the corresponding non-parametric model
%           BestIteration - Iteration number that provided best estimate
%           LastNorm - Norm of parameter change in the last iteration
%
% Example: 
%    m = wnsfBJ(iddata(y,u),[1 1 1 1],'UnstructuredOrders',[20 40 60])
%    Estimates first-order BJ model with y as input and u as output, trying
%    non-parametric ARX orders 20, 40, and 60.
%
% Other m-files required: wnsfARX.m, wnsfBJWLS.m, wnsfBJWLS_noH.m
% See also: iddata, wnsfOE

% Author: Miguel Galrinho
% Work address: Osquldas väg 10, Stockholm, Sweden
% email: galrinho@kth.se
% Created: August 2015; Last revision: 27-June-2017

% Initialization
u = z.u;
y = z.y;
N = length(y);
n_max = floor(N/4);

% Check Inputs
p = inputParser;
defaultUnstructuredOrders = round(linspace(20,n_max,10));
defaultIterations = 20;
defaultTol = 1e-4;
defaultInitialConditions = 'truncate';
expectedInitialConditions = {'zero','truncate'};
addRequired(p,'u',@(x) validateattributes(x,{'numeric'},{'2d','finite','real'}));
addRequired(p,'y',@(x) validateattributes(x,{'numeric'},{'2d','finite','real'}));
addRequired(p,'m',@(x) validateattributes(x,{'numeric'},{'vector','integer','nonnegative'}));
addParameter(p,'UnstructuredOrders',defaultUnstructuredOrders,...
    @(x) validateattributes(x,{'numeric'},{'vector','integer','positive'}));
addParameter(p,'Iterations',defaultIterations,...
    @(x) validateattributes(x,{'numeric'},{'scalar','integer','nonnegative'}));
addParameter(p,'Tol',defaultTol,...
    @(x) validateattributes(x,{'numeric'},{'scalar','positive'}));
addParameter(p,'InitialConditions',defaultInitialConditions,...
    @(x) any(validatestring(x,expectedInitialConditions)));
parse(p,u,y,m,varargin{:});

% Initializations 2
n_grid = sort(p.Results.UnstructuredOrders,'descend');
it = p.Results.Iterations;
Tol = p.Results.Tol;
InitialConditions = p.Results.InitialConditions;
models = cell(length(n_grid),3);

if strcmpi(InitialConditions,'truncate')==1
    n_max = n_grid(1);
else
    n_max = 0;
end
mf = m(1);
ml = m(2);
if length(m)==4
    mc = m(3);
    md = m(4);
    for k=1:length(n_grid)
        n = n_grid(k);
        [eta,R,L] = wnsfARX(u,y,n,InitialConditions);
        [G,H,V,it_end,it_best,norm_end] = wnsfBJWLS(eta,R,L,mf,ml,mc,md,it,Tol,u,y,n_max);
        models{k,1} = G;
        models{k,2} = H;
        models{k,3} = V;
        models{k,4} = [it_end it_best norm_end];
    end
else
    for k=1:length(n_grid)
        n = n_grid(k);
        [eta,R] = wnsfARX(u,y,n,InitialConditions);
        if k==1
            H = tf(1,[1 eta(1:n)'],1,'variable','z^-1');
        end
        [G,V,it_end,it_best,norm_end] = wnsfBJWLS_noH(eta,R,mf,ml,it,Tol,u,y,H,n_max);
        models{k,1} = G;
        models{k,2} = H;
        models{k,3} = V;
        models{k,4} = [it_end it_best norm_end];
    end
end

Vall = cell2mat(models(:,3));
[Vmin,Vmin_i] = min(Vall);
my_model.G = models{Vmin_i,1};
my_model.H = models{Vmin_i,2};
my_model.CostFunction = Vmin;
my_model.BestUnstructuredOrder = n_grid(Vmin_i);
my_model.LastIteration = models{Vmin_i,4}(1);
my_model.BestIteration = models{Vmin_i,4}(2);
my_model.LastNorm = models{Vmin_i,4}(3);

default_model = idpoly(1,my_model.G.num{1},my_model.H.num{1},my_model.H.den{1},my_model.G.den{1});

end