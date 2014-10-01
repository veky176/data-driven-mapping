% generate a SPD matrix
% A = generateSPDmatrix(1000, 1, 2);
%
% to save : dlmwrite ('input_1000.txt', A, ' '); for n = 1000
% to save : dlmwrite ('input2_1000.txt', A, ' '); for n = 1000
% dlmwrite ('input_4.txt', A, 'delimiter', ' ', 'precision', 10); for n = 4
%
function A = generateSPDmatrix(n)
% Generate a dense n x n symmetric, positive definite matrix

Q = randn(n,n);
eigen_mean = 2;
A = Q' * diag(abs(eigen_mean+randn(n,1))) * Q;
A = A + A';

L = chol (A)';
invL = inv (L);
invA = inv (A);

fprintf ('det of A = %f\nA in [%f, %f]\n', det(A), max(max(A)), min(min(A)));
fprintf ('det of L = %f\nL in [%f, %f]\n', det(L), max(max(L)), min(min(L)));
fprintf ('det of invL = %f\ninvL in [%f, %f]\n', det(invL), max(max(invL)), min(min(invL)));
fprintf ('det of invA = %f\ninvA in [%f, %f]\n', det(invA), max(max(invA)), min(min(invA)));

return;



function A = method_1 (n, min, max)
% A = rand (n, n);
A = min + (max - min) .* rand(n,n); % generate a random n x n matrix

% construct a symmetric matrix using either
% A = A + A';
A = A*A';
% The first is significantly faster: O(n^2) compared to O(n^3)

% since A(i,j) < 1 by construction and a symmetric diagonally dominant matrix
%   is symmetric positive definite, which can be ensured by adding nI

A = A + n * eye (n);
return;

function A = method_2 (n, min, max)
A = min + (max - min) .* rand(n,n); % generate a random n x n matrix
A = A*A';

number = 1;

rs = sum (A, 1)';
cs = sum (A, 2);

for i = 1:n
    value1 = rs (i);
    value2 = cs (i);
    value1 = ceil (abs (value1));
    value2 = ceil (abs (value2));
    if (value1 > value2)
        value = value1;
    else
        value = value2;
    end
    if (value > number)
        number = value;
    end
end

number = number + 1;
display (number);

A = A + number * eye(n);
return;