%% Compute clusters from Adjacency Matrix.
% Downloaded from:

% Candelier R. Adjacency matrix to clusters
%(https://la.mathworks.com/matlabcentral/fileexchange/60676-
%adjacency-matrix-to-clusters), MATLAB Central File
%Exchange. Accessed november 7, 2025.


function C = adj2cluster(A)

% Symmetrize adjacency matrix
S = A + A';

% Reverse Cuthill-McKee ordering
r = fliplr(symrcm(S));

% Get the clusters
C = {r(1)};
for i = 2:numel(r)
    if any(S(C{end}, r(i)))
        C{end}(end+1) = r(i);
    else
        C{end+1} = r(i);
    end
end