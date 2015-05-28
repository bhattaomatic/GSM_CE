function output = toepmat(input,L)
% m = 1:1:26;
% L = 5;
N = length(input);

M = zeros(N+L-1,L);

for j = 1:L
    for i = j+1:N+L
        if i-j > N
            j=j+1;
        end
        M(i-1,j) = input(i-j);
    end
end
output = M;