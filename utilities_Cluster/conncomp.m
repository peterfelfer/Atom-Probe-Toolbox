function [S,C] = conncomp(G)
  % CONNCOMP Drop in replacement for graphconncomp.m from the bioinformatics
  % toobox. G is an n by n adjacency matrix, then this identifies the S
  % connected components C. This is also an order of magnitude faster.
  %
  % [S,C] = conncomp(G)
  %
  % Inputs:
  %   G  n by n adjacency matrix
  % Outputs:
  %   S  scalar number of connected components
  %   C 
  % This function is thankfully copied from
  % https://www.alecjacobson.com/weblog/?p=4203 on the 19.03.2021
  % the code is part of the gptoolbox and written by Alex Jacobson
  [p,q,r] = dmperm(G+speye(size(G)));
  S = numel(r)-1;
  C = cumsum(full(sparse(1,r(1:end-1),1,1,size(G,1))));
  C(p) = C;
end
