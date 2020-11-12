function [x, resnorm, residual] = levmarq(func, x0)
%Parameters
TOL = 1e-3; %tolerence
u = 100; %the damping parameter mu
maxIter = 1000; %max nbr of iterations
x = x0; %initial guess vector
if(size(x, 2) > size(x, 1))
    x = x;  %must be given as column vector
end

N = length(x);
I = eye(N, N);

%counters
iter = 0;
k = 1;
xIterations(:, k) = x;

%trust region parameters
epsilon2 = 3/4;
epsilon1 = 1/4;

%choose Jacobian from user if '1' else calculate Jacobian
jacobian = 0;

if(jacobian == 1) %user gives the Jacobian
  [residual, grad] = feval(func, x);
  if(size(grad, 1) == N)
      grad = grad;
  end
else
  % else do numerical calculation of Jacobian
  residual = feval(func, x);
  grad = calcJac(func, residual, x);
end
%
scale
f(k) = (residual*residual)/2;
norm cond(k) = norm(grad*residual);

%Main loop, loop until the tolerence is reached (as long as the nbr of
%iterations are less than the maximum nbr of iterations).
while((norm cond(k) > TOL) && (iter < maxIter))
  %time-stepping parameters
  iter = iter + 1;
  k = k + 1;
  x k = x;
  %Get p using the Levenberg-Marquardt algorithm
  p = ((grad'*grad) + (u*I))\(-grad'*residual)
  x = x k + p;%update x
  if(jacobian == 1)
  [residual, grad] = feval(func, x);
  f(k) = (residual*residual)/2;
  else
  residual = feval(func, x);
  grad = calcJac(func, residual, x);
  end
  %check if the solution is good
  %calculate the ratio of the actual and predicted reductions in error:
  rho=(residual*residual-(feval(func,(x-p))*feval(func,(x-p))))/...
      ((residual*grad*p)-(p*(grad*grad+(u*I))*p));

  f(k)=(residual*residual)/2;

  % update mu with regard to the trust region
  if abs(rho) >= epsilon2
      u=u/9; %not enough damping expand the trust region (reduce mu)
  elseif abs(rho) <= epsilon1
      x=x k;
      if(jacobian==1)
          [residual,grad]=feval(func, x);
      else
          residual=feval(func, x);
          grad=calcJac(func, residual, x);
      end
      f(k)=(residual*residual)/2;
      u=11*u; %shrink the trust region (increase mu)
  end
  norm cond(k)=norm(grad*residual);
  xIterations(:,k)=x;
end
%Display nbr of iterations needed for main loop.
%disp('Iterations:')
%disp(iter)

x = xIterations(:, end);
resnorm = 2*f(end);

end

function [r, grad] = residualfunc(x)
y = [6.8 3.0 1.5 0.75 0.48 0.25 0.2 0.15]';
t = [0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0]';
r = x(1)*exp(x(2)*t)-y;

%supply gradient to be used
if(nargout > 1)
    grad = [exp(x(2)*t) x(1)*t.*exp(x(2)*t)];end
end

function grad = calcJac(func, r, x)
    n = length(x);
    grad = zeros(length(r), n);

    %The Jacobian is approximated using (backwards) finite differences
    for i = 1:n
        dx=x(i)*1e-5;
        x1=x;
        x1(i)=x1(i)+dx;
        f x = feval(func,x1);
        grad(:,i) = (f x-r)/dx;
    end
end

%
% Generate CIR distribuion from noncentral chi-square distribution.
%
% dX = a(b - X)dt + sigma*sqrt(X) dW
%
% The CIR pdf is related to the Chi distribution with 2q+2 degrees
% of freedom and non-centrality parameter
% Usage:
% Density = cirpdf(X t , t, X {t-1},t-1, a, b, sigma);
% where (X t, t) is the evaluation point at the final time
% given that the process started in (X {t-1}, t-1)
%
function res = cirpdf(Params, Model)
Data = Model.Data;
y = Data(2:end);
x = Data(1:end-1);
t = Model.time;
a = (Params(1));
b = (Params(2));
sigma = (Params(3));

d = exp(-a*(t));
c = 2*a/(sigma.ˆ2.*(1-d));
q = 2*a*b/(sigma.ˆ2)-1;

% Transformed variable:
z = 2*c*y;

% Non centrality:
lambda = 2*c*x*d;

% Degrees of freedom:
dg = 2*q+2;

% Check feller condition
%feller = 2*a*b-sigma.ˆ2;
%if(feller<0)
if(not(2*a*b>sigma.ˆ2))
      % warning('Feller condition not fullfilled')
res = 1e-100; % Return a very small probability
else
  pdf = 2*c*ncx2pdf(z,dg,lambda);
  res = sum(-log(pdf)); %return log-liklihood function

end
