% Title:            Tutorial for use of Bayesian Files
% Author:           Alex Josephson
% Last Modified:    October 2016

% This is a brief tutorial for how the Bayesian files in this directory are
% used for model calibration. This material is presented as supplemental
% material for the article: Modeling Soot Consumption with Bayesian
% Statistics.

% In this brief tutorial we will calibrate a quadratic function to a set of
% synthetic data produced.

%% Create our synthetic data

% Original model ('The correct answer') y = 2x^2 + 3x + 15
% All synthetic data are Gaussianly distributed, with a standard deviation
% of 100, along the 'correct answer'.

True_Parameters= [2 3 15];
x = 20*rand(20,1);
y = normrnd(quadratic_model(x,True_Parameters),100);

figure,
plot(x,y,'.')

%% Define our model

% Normal quadratic model with parameters A, B, and C

% model defined in the model.m file

%% Prepare all available information for evaluation by Bayesian function

% First we establish a possible range over A, B, and C that we want test.
% In this case we will test A from -5 to 5, B from -15 to 15, and C from
% -100 to 100.

Parameter_Space=[   linspace(-5,5,50);
                    linspace(-15,15,50);
                    linspace(-100,100,50)];

% Next, we establish an uninformed prior, meaning all possible combinations
% of parameters are equally likely since we have no previous notion what A,
% B, or C should be.

Prior = 1/50^3*ones(50,50,50);

% Lastly, before we run the Bayesian function we need to combine our data
% in a way that it can be read.

Data = [ y, x];

%% Run Bayesian function to get our probability space

% This could take awhile depending on the speed of your computer

[Posterior] = bayesian(Parameter_Space,Prior,Data,@quadratic_model);

%% Interpet the probability space

% Marginalize the probability space to 3 1-D arrays to form a PDF for
% each parameter.

Marginal_Space=zeros(3,size(Parameter_Space,2));
for i=1:size(Posterior,2)
    for j=1:size(Posterior,2)
        for k=1:size(Posterior,2)
            Marginal_Space(1,i)=Marginal_Space(1,i)+Posterior(i,j,k);
            Marginal_Space(2,i)=Marginal_Space(2,i)+Posterior(j,i,k);
            Marginal_Space(3,i)=Marginal_Space(3,i)+Posterior(j,k,i);
        end
    end
end
sum=zeros(3,1);
for i=1:3
    for j=1:size(Marginal_Space,2)
        sum(i)=sum(i)+Marginal_Space(i,j);
    end
end
for i=1:3
    for j=1:size(Marginal_Space,2)
        Marginal_Space(i,j)=Marginal_Space(i,j)/sum(i);
    end
end

% We can also form contours which show how the parameters are correlated to
% each other

contour_1=zeros(50,50);
contour_2=zeros(50,50);
contour_3=zeros(50,50);
for i=1:50
    for j=1:50
        for k=1:50
            contour_1(i,j)=contour_1(i,j)+Posterior(i,j,k);
            contour_2(i,j)=contour_2(i,j)+Posterior(i,k,j);
            contour_3(i,j)=contour_3(i,j)+Posterior(k,i,j);
        end
    end
end

% Visualize the resultant PDFs and contours

figure,
subplot(3,3,1)
plot(Parameter_Space(1,:),Marginal_Space(1,:),'linewidth',3)
ylabel('Marginal Posterior','FontSize',25)
xlabel('A','FontSize',25)
subplot(3,3,4)
contour(Parameter_Space(1,:),Parameter_Space(2,:),contour_1')
ylabel('B','FontSize',25)
xlabel('A','FontSize',25)
subplot(3,3,5)
plot(Parameter_Space(2,:),Marginal_Space(2,:),'linewidth',3)
ylabel('Marginal Posterior','FontSize',25)
xlabel('B','FontSize',25)
subplot(3,3,7)
contour(Parameter_Space(1,:),Parameter_Space(3,:),contour_2')
xlabel('A','FontSize',25)
ylabel('C','FontSize',25)
subplot(3,3,8)
contour(Parameter_Space(2,:),Parameter_Space(3,:),contour_3')
xlabel('B','FontSize',25)
ylabel('C','FontSize',25)
subplot(3,3,9)
plot(Parameter_Space(3,:),Marginal_Space(3,:),'linewidth',3)
ylabel('Marginal Posterior','FontSize',25)
xlabel('C','FontSize',25)

%% Validate

% Take the mode of each PDF to determine the optimal parameter value to be
% used in our model.


parameters=zeros(1,3);

% Index the maximum value
parametersi=ones(1,3);
for i=1:size(Marginal_Space,2)
    if Marginal_Space(1,i)>Marginal_Space(1,parametersi(1));
        parametersi(1)=i;
    end
    if Marginal_Space(2,i)>Marginal_Space(2,parametersi(2));
        parametersi(2)=i;
    end
    if Marginal_Space(3,i)>Marginal_Space(3,parametersi(3));
        parametersi(3)=i;
    end
end

% Estimate maximum value using a quadratic regression
A=zeros(3,3,3);
B=zeros(3,1,3);
Sol=zeros(3,3);
for i=1:3
    for j=1:3
        temp=-2+j;
        A(j,1,i)=(Parameter_Space(i,parametersi(i)+temp))^2;
        A(j,2,i)=Parameter_Space(i,parametersi(i)+temp);
        A(j,3,i)=1;
        B(j,1,i)=Marginal_Space(i,parametersi(i)+temp);
    end
    Sol(i,:)=A(:,:,i)\B(:,:,i);
end

% Take a derivative and solve for maximum value
for i=1:3
    parameters(i)=-Sol(i,2)/(2*Sol(i,1));
end

% Check how our calibration did against the 'correct answer'

x_modeled = linspace(0,20,100);
y_modeled = parameters(1)*x_modeled.^2+parameters(2)*x_modeled+parameters(3);

figure,
plot(x,y,'.',x_modeled,y_modeled)

disp('True Parameters')
disp(True_Parameters)
disp('Calibrated Parameters')
disp(parameters)

%% Add more data

% Let's say we run accross more data that we want to include in our model
% calibration. It's really simple, we simply run our Bayesian function
% again using the previous posterior as our new prior along with the new
% data.

% New synthetic data

x_new = 50+25*rand(20,1);
y_new = normrnd(quadratic_model(x_new,True_Parameters),500);

figure,
plot(x,y,'.',x_new,y_new,'.')

%% Run Bayesian function again using previous Posterior as our new Prior

Data_new = [ y_new, x_new];

% This could take awhile depending on the speed of your computer. The first
% expression has no significant value but it keeps the likelihood values in
% the bayesian function large enough to avoid numerical truncation errors.
Posterior = 10^200 * Posterior;

[Posterior] = bayesian(Parameter_Space,Posterior,Data_new,@quadratic_model);

%% Then we visualize everything again

Marginal_Space=zeros(3,size(Parameter_Space,2));
for i=1:size(Posterior,2)
    for j=1:size(Posterior,2)
        for k=1:size(Posterior,2)
            Marginal_Space(1,i)=Marginal_Space(1,i)+Posterior(i,j,k);
            Marginal_Space(2,i)=Marginal_Space(2,i)+Posterior(j,i,k);
            Marginal_Space(3,i)=Marginal_Space(3,i)+Posterior(j,k,i);
        end
    end
end
sum=zeros(3,1);
for i=1:3
    for j=1:size(Marginal_Space,2)
        sum(i)=sum(i)+Marginal_Space(i,j);
    end
end
for i=1:3
    for j=1:size(Marginal_Space,2)
        Marginal_Space(i,j)=Marginal_Space(i,j)/sum(i);
    end
end

contour_1=zeros(50,50);
contour_2=zeros(50,50);
contour_3=zeros(50,50);
for i=1:50
    for j=1:50
        for k=1:50
            contour_1(i,j)=contour_1(i,j)+Posterior(i,j,k);
            contour_2(i,j)=contour_2(i,j)+Posterior(i,k,j);
            contour_3(i,j)=contour_3(i,j)+Posterior(k,i,j);
        end
    end
end

% Visualize the resultant PDFs and contours

figure,
subplot(3,3,1)
plot(Parameter_Space(1,:),Marginal_Space(1,:),'linewidth',3)
ylabel('Marginal Posterior')
xlabel('A','FontSize',25)
subplot(3,3,4)
contour(Parameter_Space(1,:),Parameter_Space(2,:),contour_1')
ylabel('B','FontSize',25)
xlabel('A','FontSize',25)
subplot(3,3,5)
plot(Parameter_Space(2,:),Marginal_Space(2,:),'linewidth',3)
ylabel('Marginal Posterior','FontSize',25)
xlabel('B','FontSize',25)
subplot(3,3,7)
contour(Parameter_Space(1,:),Parameter_Space(3,:),contour_2')
xlabel('A','FontSize',25)
ylabel('C','FontSize',25)
subplot(3,3,8)
contour(Parameter_Space(2,:),Parameter_Space(3,:),contour_3')
xlabel('B','FontSize',25)
ylabel('C','FontSize',25)
subplot(3,3,9)
plot(Parameter_Space(3,:),Marginal_Space(3,:),'linewidth',3)
ylabel('Marginal Posterior','FontSize',25)
xlabel('C','FontSize',25)

for i=1:size(Marginal_Space,2)
    if Marginal_Space(1,i)>Marginal_Space(1,parametersi(1));
        parametersi(1)=i;
    end
    if Marginal_Space(2,i)>Marginal_Space(2,parametersi(2));
        parametersi(2)=i;
    end
    if Marginal_Space(3,i)>Marginal_Space(3,parametersi(3));
        parametersi(3)=i;
    end
end

for i=1:3
    for j=1:3
        temp=-2+j;
        A(j,1,i)=(Parameter_Space(i,parametersi(i)+temp))^2;
        A(j,2,i)=Parameter_Space(i,parametersi(i)+temp);
        A(j,3,i)=1;
        B(j,1,i)=Marginal_Space(i,parametersi(i)+temp);
    end
    Sol(i,:)=A(:,:,i)\B(:,:,i);
end

for i=1:3
    parameters(i)=-Sol(i,2)/(2*Sol(i,1));
end

% Check how our calibration did against the 'correct answer'

x_modeled = linspace(0,80,80);
y_modeled = parameters(1)*x_modeled.^2+parameters(2)*x_modeled+parameters(3);

figure,
plot(x,y,'.',x_new,y_new,'.',x_modeled,y_modeled)

disp('True Parameters')
disp(True_Parameters)
disp('Calibrated Parameters')
disp(parameters)

% Was there improvement with the additional data?
% Any time more data is added, the model will get better, eventually, with
% enough data we would recover the 'correct answer'.