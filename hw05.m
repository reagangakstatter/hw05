% Author: Reagan Gakstatter / reg0052@auburn.edu
% Date: 2024-11-08
% Assignment Name: hw05

classdef hw05
    methods (Static)
        
        function ret = p1(data, powers)
            % Implementation of Richardson Extrapolation
            % Assume the expansion of f(h) = f(0) + c_1 h^{alpha_1} + c_2 h^{alpha_2} + ... + c_n h^{alpha_n} + ...
            %
            %:param: data: a vector of values f(2^(-i) h), i=1,2,...,n
            %:param: powers: a vector of powers (alpha_1, ..., alpha_{n-1})
            %
            %:return: the extrapolated value of f(0) using Richardson Extrapolation
            
            m = length(data);
            R = zeros(m, m);  % Initialize Richardson table
            R(:, 1) = data(:); % First column is initial data values
            
            % Write your code here.
            for j = 2:m
                for i = 1:(m - j + 1)
                    exponent = powers(j - 1);
                    R(i, j) = (2^exponent * R(i+1, j-1) - R(i, j-1)) / (2^exponent - 1);
                end
            end
            
            % Best approximation of f(0) is in the first element of the last column
            ret = R(1, m);
        end

        function R = p2(beta)
            % Compute the value of the series 
            %   sum_{k=0}^(\infty) ((-1)^k /(2k + 1)^{beta})

            %:param: beta: a real number on (0, 1].
            %:return: R: the value of the series
            
            m = 15;  % Use m = 15 to limit rounding errors
            data = zeros(1, m);

            % Write your code here.
            % Calculate f(h, beta) for each term using different step sizes h = 2^(-k)
            for k = 1:m
                h = 2^(-k);
                f_h_beta = 0;
                % Compute the partial sum for the series up to 1/h terms
                for i = 0:1/h
                    f_h_beta = f_h_beta + ((-1)^i) / (2*i + 1)^beta;
                end
                data(k) = f_h_beta;
            end
            
            % Vector of powers for Richardson extrapolation
            powers = beta + (0:m-2);
            
            % Apply Richardson extrapolation to get the best approximation
            R = hw05.p1(data, powers);
        end

        function coefs = p3(shifts)
            % Compute the coefficients of the finite difference scheme for f'(x)
            % using the formula
            % f'(x) \approx \frac{1}{h} (c_0 f(x_0) + c_1 f(x_1) + ... + c_n f(x_n)) + O(h^n)

            %:param: shifts: a vector of shifts (a_0, a_1, ..., a_n), the nodes are x_i = x + a_i h
            %:return: coefs: a vector of coefficients (c_0, c_1, ..., c_n)

            m = length(shifts);
            coefs = zeros(m, 1);
            
            % Compute Lagrange interpolation derivative at each shift to get coefficients
            for k = 1:m
                L_prime = 0;
                for n = 1:m
                    if n ~= k
                        prod_term = 1 / (shifts(k) - shifts(n));
                        for j = 1:m
                            if j ~= k && j ~= n
                                prod_term = prod_term * (0 - shifts(j)) / (shifts(k) - shifts(j));
                            end
                        end
                        L_prime = L_prime + prod_term;
                    end
                end
                coefs(k) = L_prime;
            end
        end
        
    end
end