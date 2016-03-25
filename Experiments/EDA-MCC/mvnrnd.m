function r = mvnrnd(mu,sigma,cases);
%MVNRND Random matrices from the multivariate normal distribution.
%   R = MVNRND(MU,SIGMA,CASES) returns a matrix of random numbers chosen   
%   from the multivariate normal distribution with mean vector, MU, and 
%   covariance matrix, SIGMA. CASES is the number of rows in R.
%
%   SIGMA is a symmetric positive definite matrix with size equal to the 
%   length of MU or simply a scalar.

% Scott J Gaffney   5 February 2003
% Department of Information and Computer Science
% University of California, Irvine.

PROGNAME = 'mvnrnd';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end


if (~isvector(mu))
   error('MU must be a vector.');
end
mu = mu(:);
d = length(mu);
if (exist('cases')~=1 | isempty(cases))
  cases = 1;
end

[m,n] = size(sigma);
if (~(m*n==1) & (m~=n | m~=d))
   error('SIGMA must be square and match the length of MU.');
end

[T,p] = chol(sigma);
if (p~=0)
  r = nan;
  return;
end

mu = mu(:,ones(cases,1));  % repmat
r = randn(cases,d)*T + mu';

