%% Subroutine to compute path length
%% Source Unknown!!

function [L,Deg,Nedges] = characteristic_path_length(A)
    % Convert to graph object and use built-in distance calculator
    G = graph(A);
    Deg=degree(G);
    Nedges=numedges(G);
    D = distances(G);  % All-pairs shortest paths
    D(D == inf) = 0;   % Disconnected components -> 0 (exclude from mean)
    
    n = size(D, 1);
    valid_distances = D(triu(true(n), 1));  % Upper triangle, exclude diagonal
    valid_distances = valid_distances(valid_distances > 0);  % Remove zeros (disconnected/inf)
    
    L = mean(valid_distances);
end