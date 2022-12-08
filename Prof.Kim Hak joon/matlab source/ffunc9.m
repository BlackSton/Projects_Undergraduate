%% Asymmetric Pseudo Voigt Function (1st Derivative)
function y = ffunc9(x,varargin)
%% Read Function Parameters
% This function uses a variable amount of parameters, considering that
% every additional peak will add a set of five paramerts to the function.
varargin = cell2mat(varargin);
c = varargin(1:5:end); % peak centers
a = varargin(2:5:end); % peak heights
w = varargin(3:5:end); % peak widths
f = varargin(4:5:end); % shape factor
s = varargin(5:5:end); % asymmetry factor
x = repmat(x,1,size(c,2)); % wavenumbers

%% Substitute Parameters
% In this step, some calculations that are performed multiple times are out
% sourced in order to decrease computing operations.
x = x-c;
%    [~,idx] = max(x==0,[],1);
%    maxX = w*30;
%    x(x>maxX) = NaN;
%    x(x<-maxX) = NaN;
u = exp(s.*x);
v = u.*x;
z = w.^2;
t = u+s.*v+1;

%% Function
% y = -nansum( ...
y = -sum( ...
    + (24*a.*f.*z.*x.*(u + 1).*t)./(v.^2 + 2.*v.*x + 12.*z + x.^2).^2 ... Lorentzian 
    - (a.*x.*exp(-(x.^2.*(u + 1).^2)./(8.*z)).*(f - 1).*(u + 1).*t)./(4.*z), ... Gaussian
    2);
end % closing function