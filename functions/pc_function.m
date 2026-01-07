%Function for estimating p_c, as described in 
% Radicci F, "Predicting percolation thresholds in networks", Phys. Rev. E 91, 010801(R) (2015)

function [p_c_estimates, results] = pc_function(A, varargin)
% COMPUTE_PERCOLATION_THRESHOLD Estimate percolation threshold using multiple methods

% Input:
%   A - adjacency matrix (symmetric, undirected graph)
%   varargin - optional parameters

% Output:}

%   p_c_estimates - struct containing estimates from different methods
%   results - detailed results from each method

% Methods:
%   A) Inverse of largest eigenvalue
%   B) Largest eigenvalue of non-backtracking matrix
%   C) First and second moments of degree distribution


p = inputParser;

addParameter(p, 'verbose', true, @islogical);
addParameter(p, 'nonbacktracking', false, @islogical);   % <-- agregado

parse(p, varargin{:});

verbose = p.Results.verbose;
nonbacktracking = p.Results.nonbacktracking;



% Validate input

if ~issymmetric(A)
    
    warning('Adjacency matrix is not symmetric. Making it symmetric.');
    
    A = max(A, A');
    
end



if any(diag(A) ~= 0) % remove self-loops
    
    A = A - diag(diag(A));
    
end


A = sparse(A);
n = size(A, 1);

% Initialize results structure
p_c_estimates = struct();
results = struct();

%% Method A: Inverse of largest eigenvalue

% Compute largest eigenvalue

opts.tol = 1e-8;
opts.maxit = 1000;
lambda_max = eigs(A, 1, 'la', opts); % 'la' = largest algebraic

p_c_A = 1 / lambda_max;
results.method_A.lambda_max = lambda_max;
results.method_A.p_c = p_c_A;
p_c_estimates.spectral = p_c_A;


%% Method B: Non-backtracking matrix
if nonbacktracking
    
    try
        
        % Build non-backtracking matrix
        
        B = build_non_backtracking_matrix(A);
        
        % Compute largest eigenvalue
        
        if size(B, 1) > 0
            
            lambda_nb = eigs(B, 1, 'la', opts);
            
            p_c_B = 1 / lambda_nb;
            
        else
            
            lambda_nb = 0;
            
            p_c_B = Inf;
            
        end
        
        results.method_B.non_backtracking_matrix = B;
        results.method_B.lambda_nb = lambda_nb;
        results.method_B.p_c = p_c_B;
        
        p_c_estimates.non_backtracking = p_c_B;
        
        
        
    catch ME
        
        warning('Non-backtracking method failed: %s', ME.message);
        
        results.method_B.error = ME.message;
        
        p_c_estimates.non_backtracking = NaN;
        
    end
    
end

%% Method C: Degree Distribution Moments

% Compute degrees

degrees = full(sum(A, 2));

% Remove isolated nodes for moment calculation

valid_degrees = degrees(degrees > 0);



if length(valid_degrees) >= 2
    
    % First moment (mean degree)
    
    k_mean = mean(valid_degrees);
    % Second moment
    
    k2_mean = mean(valid_degrees.^2);
    % Estimate using Molloy-Reed criterion approximation
    
    % p_c â‰ˆ <k> / (<k^2> - <k>)
    
    if k2_mean > k_mean
        
        p_c_C = k_mean / (k2_mean - k_mean);
        
    else
        p_c_C = Inf;
    end
    % Alternative: p_c â‰ˆ 1 / (k2_mean/k_mean - 1)
    p_c_C_alt = 1 / (k2_mean/k_mean - 1);
    
    results.method_C.degrees = degrees;
    results.method_C.k_mean = k_mean;
    results.method_C.k2_mean = k2_mean;
    results.method_C.p_c_molloy_reed = p_c_C;
    results.method_C.p_c_alternative = p_c_C_alt;
    
    p_c_estimates.degree_moments = p_c_C;
    p_c_estimates.degree_moments_alt = p_c_C_alt;
    
else
    warning('Not enough degree information for moment method');
    results.method_C.error = 'Insufficient degree distribution';
    p_c_estimates.degree_moments = NaN;
    p_c_estimates.degree_moments_alt = NaN;
end
end
