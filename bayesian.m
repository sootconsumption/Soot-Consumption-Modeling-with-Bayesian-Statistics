function [posterior] = bayesian(sample_space,prior,data,model)

% Note: The inputs to this function should include a sample space, a
% prior,the model to be calibrated and the data over which the model should
% be tested. The sample space should be a two dimensional matrix with the
% rows representing different parameters and the columns representing the
% different values over which those parameters will be tested. The prior
% should be a multi-dimensional matrix with the number of dimensions equal
% to the number of parameters to be tested. Data should be a two
% dimensional matrix containing all the data to be in comparison in the
% first column and all dependent data to be used in the model in the other
% columns.

n = size(sample_space,2);
p = size(sample_space,1);

% Establish Jeffrey's Prior for the prior of the Sigma

sigma = linspace(1,100,n);
sigma_prior = zeros(n,1);
for i = 1:n
    sigma_prior(i) = (2/sigma(n)^2)^0.5;
end

% Discretize the data matrix into dependent to be compared and independent
% data which are used in the model.

dependent_data = data(:,1);
data(:,1) = [];

% Initiate matrices

posterior = zeros(size(prior));
parameters = zeros(p,1);
sum = 0;

% Perform the Bayesian Process

for i = 1:n^p
    temp = i;

    % Determine what combination of parameters to be tested
    for j = 1:p
        k = p-j+1;
        parameters(k) = sample_space(k,ceil(temp/n^(k-1)));
        if temp > n^(k-1)
             temp = temp-floor(temp/n^(k-1))*n^(k-1);
             if temp == 0
                 temp = 0.1;
             end
        end
    end
    
    % Vary over different values of sigma
    
    for k = 1:n
        predicted = model(data,parameters);
        likelihood = prod((exp(-(dependent_data(:)-predicted).^2/(2*sigma(k)^2))/(sigma(k)*(2*pi)^0.5)));
        posterior(i) = posterior(i)+prior(i)*sigma(k)*likelihood; 
    end
    sum = sum+posterior(i);
end

% Normalize the PDF
for i = n^p
    posterior(i) = posterior(i)/sum;
end
end