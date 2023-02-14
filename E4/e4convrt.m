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

function e4convrt(file)

NV = [...
'EEOPTION';...
'ee_dvee ';...
'ee_dv   ';...
'ee2thd  ';...
'eedvper ';...
'eeinit  ';...
'eelnsrch';...
'eemin   ';...
'estr_dv ';...
'estr2thd';...
'fvfast  ';...
'fvgarch ';...
'fvmiss  ';...
'fvmod   ';...
'fvper   ';...
'fvsc    ';...
'garch2ee';...
'sc2ee   ';...
'seteeopt';...
'sifmiss ';...
'sifmod  ';...
'sifres  ';...
'sifriss ';...
'sifrsc  ';...
'sifsc   ';...
'thd2eee ';...
'thd2ee  ';...
'thd2estr';...
'vecee   '];
NN = [...
'E4OPTION';...
'ss_dvss ';...
'ss_dv   ';...
'ss2thd  ';...
'ssdvper ';...
'e4init  ';...
'e4lnsrch';...
'e4min   ';...
'str_dv  ';...
'str2thd ';...
'lffast  ';...
'lfgarch ';...
'lfmiss  ';...
'lfmod   ';...
'lfper   ';...
'lfsc    ';...
'garch2ss';...
'sc2ss   ';...
'sete4opt';...
'fismiss ';...
'fismod  ';...
'fisres  ';...
'fisriss ';...
'fisrsc  ';...
'fissc   ';...
'thd2ssss';...
'thd2ss  ';...
'thd2str ';...
'vecss   '];

N = 0;
for i=1:size(NV,1)
    N = N + replace(file,file,deblank(NV(i,:)),deblank(NN(i,:)));
    N = N + replace(file,file,deblank(upper(NV(i,:))),deblank(upper(NN(i,:))));
end

fprintf(1,'N�mero total de cambios: %3d\n',N);