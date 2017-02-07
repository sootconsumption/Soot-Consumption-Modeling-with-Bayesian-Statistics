function [y] = quadratic_model(x,parameter)
y = parameter(1)*x.^2+parameter(2)*x+parameter(3);
end