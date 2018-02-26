function [default_model,my_model] = wnsfOE(z,m,varargin)
%wnsfOE - Estimate a SISO output-error model with WNSF
%y(t) = B(q)/F(q) u(t) + e(t),
%where B(q) =     b_1 q^{-1} + ... + b_{mb} q^{-mb} 
%and   F(q) = 1 + f_1 q^{-1} + ... + f_{mf} q^{-mf}
%
% Syntax:  [model_idpoly_format,model_alternative_format] = wnsfOE(data,polynomial_orders,Options)
%
% Inputs:
%    data - Input and Output signals in iddata format
%    polynomial_orders - Vector with two entries [mf mb]
%    'UnstructuredOrders' (optional) - Followed by a vector of
%       non-parametric model orders to be tried out (default: 10 linearly
%       spaced orders between 20 and sample_size/3)
%    'Iterations' (optional) - Followed by an integer with the maximum
%       number of WLS iterations (default: 20).
%    'Tol' (optional) - Followed by float with parameter change tolerance
%       for stopping iterations (default: 1e-4)
%    'InitialConditions' (optional) - Folloed by string 'zero', 'truncate',
%       or 'estimate' on how to handle initial conditions (default:
%       'estimate')
%
% Outputs:
%    model_idpoly_format - Estimated model in idpoly format
%    model_alternative_format - Estimated model as a structure with
%    following entries:
%           G - Plant model
%           CostFunction - PEM cost function value for the estimated model
%              (estimated noise variance)
%           BestUnstructuredOrder - Non-parametric FIR model order used for
%              the final parametric estimate
%           LastIteration - Number of weighted least-squares iterations
%              performed for the corresponding non-parametric model
%           BestIteration - Iteration number that provided best estimate
%           LastNorm - Norm of parameter change in the last iteration
%
% Example: 
%    m = wnsfOE(iddata(y,u),[1 1],'UnstructuredOrders',[20 40 60])
%    Estimates first-order OE model with y as input and u as output, trying
%    non-parametric FIR orders 20, 40, and 60.
%
% Other m-files required: wnsfFIR.m, wnsfOEWLS.m
% See also: iddata, wnsfBJ

% Author: Miguel Galrinho
% Work address: Osquldas väg 10, Stockholm, Sweden
% email: galrinho@kth.se
% Created: August 2015; Last revision: 5-April-2017

% Initialization
u = z.u;
y = z.y;
N = length(y);
n_max = floor(N/3);

% Check Inputs
p = inputParser;
defaultUnstructuredOrders = round(linspace(20,n_max,10));
defaultIterations = 20;
defaultTol = 1e-4;
defaultInitialConditions = 'estimate';
defaultInitialParameters = [];
expectedInitialConditions = {'zero','truncate','estimate'};
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
addParameter(p,'InitialParameters',defaultInitialParameters,...
    @(x) validateattributes(x,{'numeric'},{'vector'}));
parse(p,u,y,m,varargin{:});

% Initializations 2
n_grid = sort(p.Results.UnstructuredOrders,'descend');
it = p.Results.Iterations;
Tol = p.Results.Tol;
theta0 = p.Results.InitialParameters;
InitialConditions = p.Results.InitialConditions;
models = cell(length(n_grid),3);

if strcmpi(InitialConditions,'truncate')==1
    n_max = n_grid(1);
elseif strcmpi(InitialConditions,'zero')==1
    n_max = 0;
else
    n_max = -1;
end
mf = m(1);
ml = m(2);
for k=1:length(n_grid)
    n = n_grid(k);
    [eta,R] = wnsfFIR(u,y,n,InitialConditions);
    [G,V,it_end,it_best,norm_end] = wnsfOEWLS(eta,R,mf,ml,it,Tol,u,y,n_max,theta0);
    models{k,1} = G;
    models{k,2} = V;
    models{k,3} = [it_end it_best norm_end];
end

Vall = cell2mat(models(:,2));
[Vmin,Vmin_i] = min(Vall);
my_model.G = models{Vmin_i,1}(1);
if n_max==-1
    my_model.Gini = models{Vmin_i,1}(2);
end
my_model.CostFunction = Vmin;
my_model.BestUnstructuredOrder = n_grid(Vmin_i);
my_model.LastIteration = models{Vmin_i,3}(1);
my_model.BestIteration = models{Vmin_i,3}(2);
my_model.LastNorm = models{Vmin_i,3}(3);

default_model = idpoly(1,my_model.G.num{1},1,1,my_model.G.den{1},Vmin);

end