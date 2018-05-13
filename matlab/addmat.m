% adding 2 matrices together
function [C, m, n] = addmat(A, B)
    % function definition signature
    % function [return_values, ...]  = func_name(args, ...)
    
    [m, n] = size(A);
    for i = 1:m
        for j = 1:n
            C(i, j) = A(i, j) + B(i, j);
        end
    end
end