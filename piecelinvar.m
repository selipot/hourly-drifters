function [v,varargout] = piecelinvar(x,y,u,varargin)
%PIECELIN  Piecewise linear interpolation.
%  v = piecelinvar(x,y,u) finds the piecewise linear L(x)
%  with L(x(j)) = y(j) and returns v(k) = L(u(k)).
%  First divided difference
% Shane Elipot, University of Miami, RSMAS, selipot@rsmas.miami.edu, January 2017

   diffx = diff(x,1,1);
   diffy = diff(y,1,1);
   delta = diffy./diffx;
%  Find subinterval indices k so that x(k) <= u < x(k+1)

   n = length(x);
   k = ones(size(u));
   for j = 2:n-1
      k(x(j) <= u) = j;
   end

%  Evaluate interpolant

   s = u - x(k);
   v = y(k) + s.*delta(k);

% Evaluate error

   if nargout == 2
       
       if nargin < 4
           error('you request output error but did you supply input errors');
       else
           deltay = varargin{1};
           varargout{1} = ((1-s./diffx(k)).^2.*deltay(k).^2 + (s./diffx(k)).^2.*deltay(k+1).^2).^0.5;
       end
   end
   