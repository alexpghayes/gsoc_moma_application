function [y] = soft_thr(x, lambda, ~)
y = sign(x) .* max(abs(x) - lambda, 0);
end
