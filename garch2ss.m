function [Phi, Gam, E, H, D, C, Q, Phig, Gamg, Eg, Hg, Dg] = garch2ss(theta, din)
% garch2ss - Converts a GARCH model in THD format to the corresponding SS form.
%    [Phi, Gam, E, H, D, C, Q, Phig, Gamg, Eg, Hg, Dg] = ...
%                                          garch2ss(theta, din)
% Only for models with GARCH structure in the error term.
%
% 7/3/97

% Copyright (C) 1997 Jaime Terceiro
% 
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2, or (at your option)
% any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details. 
% 
% You should have received a copy of the GNU General Public License
% along with this file.  If not, write to the Free Software Foundation,
% 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

[H_D, type, m, r, s, n, np, userflag, userf, isinnov, szpriv] = e4gthead(din);

if fix(type/100) ~= 1, e4error(14); end

if userflag
   [Phi, Gam, E, H, D, C, Q, Phig, Gamg, Eg, Hg, Dg] = garc2own(theta, din,userf(1,:));
   return;
end

[H_D, type, m, r, s, n, np, userflag, userf, isinnov, szpriv] = e4gthead(din(H_D+1:2*H_D));
[Phi, Gam, E, H, D, C, Q] = thd2ss(theta(1:np,:), din(H_D+1:2*H_D+szpriv(1)));
[Phig, Gamg, Eg, Hg, Dg] = thd2ss(theta(np+1:size(theta,1),:), din(2*H_D+szpriv(1)+1:size(din,1)));

