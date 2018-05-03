% function [Q,R] = MGMa(A)
% % Computation of the "skinny" QR decomposition of an m-by-n matrix
% % (m>=n) using the the Modified Gram-Schmidt algorithm.
% % See Trefethen and Bau, Algorithm 8.1
%  [m,n] = size(A);
%  if m<n, fprintf('ERROR: A should be an m-by-n matrix with m >= n.\n'); end
%  Q = zeros(m,n);
%  R = zeros(n,n);
%  Q = A;
%  for j=1:n
%     R(j,j) = norm(Q(:,j));
%     Q(:,j) = Q(:,j)/R(j,j);
%     for k=j+1:n
%         R(j,k) = Q(:,j)'*Q(:,k);
%         Q(:,k) = Q(:,k) - R(j,k)*Q(:,j);
%     end
%  end

function [Q, R] = MGMa(A)
%MGS Performs a modified Gram-Schmidt orthogonalization
%
%   [Q, R] = mgs(A)
%
% Performs a modified Gram-Schmidt orthogonalization. This is a more stable
% way to compute a QR-factorization.
%
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction v0.2b.
% The toolbox can be obtained from http://www.cs.unimaas.nl/l.vandermaaten
% You are free to use, change, or redistribute this code in any way you
% want. However, it is appreciated if you maintain the name of the original
% author.
%
% (C) Laurens van der Maaten
% Maastricht University, 2007

   
% Perform Gram-Schmidt orthogonalization
[m n] = size(A);
V = A;
R = zeros(n, n);
for i=1:n
    R(i,i) = norm(V(:,i));
    V(:,i) = V(:,i) / R(i, i);
    if (i < n)
        for j = i+1:n
            R(i,j) = V(:,i)' * V(:,j);
            V(:,j) = V(:,j) - R(i, j) * V(:,i);
        end
    end
end
Q = V;

end