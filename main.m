clc;clear all;close all;
format long
load CIRDataSet %load intrest rate data from file
Model.Data = data';
Model.t = t';

%parameters
Model.time = 2; %time step
N = length(Model.Data);
x = Model.Data(1:end-1); % Time series of observed interest rates

%Perform OLS on the discrete model to obtain initial
%parameters a,b and sigmaaaaaaaaa
y = [ ones(N-1,1) x];
fit = (y'*y)ˆ( -1)*(y'*Model.Data(2:end));
a = -log(fit(2))/Model.time;
dx = diff(Model.Data)./sqrt(x);
xx = Model.time./sqrt(x);
yy = Model.time.*sqrt(x);
vector = [xx yy];
coeff = vector\dx ;
residual = vector*coeff - dx; %get the residual
b = mean(data);
%sigma = std(residual)/Model.time;
sigma = (sqrt(norm(residual)ˆ2/(N-2)))/Model.time;

%InitialParams = [a b sigma];
%Use fminsearch to estimate parameters
%[x, resnorm, residual] = levmarq(@(Params) cirpdf(Params, Model), InitialParams)'
