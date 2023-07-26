function res = AICc_residuals( m , v , y , k  , n, type , width )% could have 6th and 7th inputs
%
% Corrected Akaike information criterion.
%
% Usage: res = AICc_residuals(m,v,y,k,n,type,width)
%
% Akaike information criterion corrected for finite sample size.
% The correction assumes that the model is univariate, linear, and
% has normally-distributed residuals.
%
% Call this function with five arguments (m,v,y,k,n) to get the normal sum-of-squares error function.
%
% INPUT:
%   m    vector of measurements
%   v    variances of the measurements
%   y    direct theory values
%   k    number of parameters
%   n    sample size
%   type - type of error-function, defaults to sum-of-squares, 't' for
%          tanh(dx^2), 'l' for log(d^2+dx^2)
%   width - width of robust cut-off, defaults to 3
% OUTPUT:
%   res - scaled-weighted residuals - to be used for Corrected Akaike
%         information criterion optimization
%
%
% IV 2016
% BG 2022
%
% Copyright I Virtanen <ilkka.i.virtanen@oulu.fi> and B Gustavsson <bjorn.gustavsson@uit.no>
% This is free software, licensed under GNU GPL version 2 or later

if nargin == 6
  width = 3;
end
if nargin > 5 && type == 't'
  % ss = sum( 2*width^2*tanh(( m - y ).^2 ./ (v*width^2*2) ),'all' ) + 2*k + 2*k*(k+1)/(n-k-1);
  % res = grad_y(ss)                                                   ^these terms constant
  res = (tanh(( m - y ).^2 ./ (v*width^2*2) ) -1).*( m - y )*2;
elseif nargin > 5 && type == 'l'
  try
    % ss = sum( 2*width^2*log(1/2*( m - y ).^2 ./ v  + width^2 ) - 2*width^2*log(width^2),'all') + 2*k + 2*k*(k+1)/(n-k-1);
    % res = grad_y(ss)                                             ^these terms constant ->
    % gss = ( 2*width^2*log(1/2*( m - y ).^2 ./ v  + width^2 ));
    res = 2*width^2*( m - y )./(v.*(width^2 + ( m - y ).^2 ./ v ) );
  catch
      keyboard
  end
else
  % ss = ( sum( ( m - y ).^2 ./ v ,'all') ) + 2*k + 2*k*(k+1)/(n-k-1);
  % res = grad_y(ss)                          ^these terms constant, have
  %                                           to add them after
  %                                           optimisation
  res = ( 2*( ( m - y )./ v ) );
end
res = res - min(y,0); % Penalty-term for negative electron-densities (IV: are these even possible?)
end