function result = test_spectral_radius(d)
% Computes the spectral radius as an approximate test of stationarity.
% This is not strictly speaking a stationarity test, but it is a good proxy
% on whether a granger causality estimation of a data matrix including this
% trial would fail.
%
% INPUT:    d should be of size  ROIS x T x 2

regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
morder=1;
[A, trash] = tsdata_to_var(d,morder,regmode);

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');
pn1 = (p-1)*n;

A1 = [reshape(A,n,p*n); eye(pn1) zeros(pn1,n)];



% calculate spectral radius
radius = max(abs(eig(A1)));

if radius >= 1
    result = 0;
else
    result = 1;
end

end

