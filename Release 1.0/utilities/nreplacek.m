function [M,m] = nreplacek(n,k)
% [M,m] =  = nreplacek(n,k)
%
% Enumerates all outcomes of drawing k times from the set [1 2 ... n] WITH
% replacement, but WITHOUT caring for order.
% m = nchoosek( n+k-1, k) is the number of outcomes
% M is a (possibly) huge m-by-k matrix, where each row is a configuration
% % EXAMPLE:
% [M,m] = nreplacek(4,2)
% 
% M =
% 
%      4     4
%      4     3
%      4     2
%      4     1
%      3     3
%      3     2
%      3     1
%      2     2
%      2     1
%      1     1
% 
% 
% m =
% 
%     10

% m = nchoosek(n+k-1,k);
% M = zeros(m,k);

% Jonathan Epperlein, July 2011
% jpe(at)engineering.ucsb.edu
% Version 0.1

if n<1||abs(round(n)-n)>eps
     error('nreplacek:iputchk',...
         'In nreplacek(n,k), n needs to be a positive integer. You supplied %f',n);
end
if k<1||k>n||abs(round(k)-k)>eps
     error('nreplacek:iputchk',...
         'In nreplacek(n,k), k needs to be a positive integer less than or equal to n. You supplied %f',k);
end

M = [];
m = 0;
for iter=n:-1:1;
    recrs(k-1, iter)
end

function recrs(count,   row)
    if count>=1
       for jter=row(end):-1:1
           recrs(count-1, [row jter]);
       end
    elseif count==0
        M = [M;row];
        m = m+1;
    else
        error('nreplacek:recrs','counter messed up, it''s value is %f',count);
    end
end
end