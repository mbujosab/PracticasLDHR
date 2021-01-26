%% usage: [AR,s2_e,n_retardos] = cib (serie,rango_retardos,IA) 
%% 
%% Calcula el orden de un polinomio autorregresivo seg�n el  
%% Criterio de Informaci�n Bayesiano 
%%  
%% Serie          > serie temporal  
%% rango_retardos > [orden m�nimo retardos, orden m�ximo] % ([1 T/3]) 
%%                    (Si escalar => ese sera el orden del modelo AR) 
%%                    (T = longitud de la muestra) 
%% IA             > Informaci�n Adicional                 % (0) 
%%                    (>0 calcula ra�ces, =2 informaci�n por pantalla)  
%% 
%% AR             > polinomio autorregresivo estimado 
%% s2_e           > varianza residual del modelo AR 
%% n_retardos     > numero de retardos empleados en el modelo AR 
%% 
%% Copyright (C) 1999, 2000, 2001 Marcos Bujosa. 
 
%% Author: Marcos Bujosa <marcos.bujosa@ccee.ucm.es> 
%% References: Box, G. E. P., Jenkins, G. M., and Reinsel, G. C. (1994).
%% Time series analysis: forecasting and control. Prentice Hall, Inc.,
 
function [AR,s2_e,n_retardos] = cib (y,rango_retardos,IA) 

  T=length(y); y=y(:); 
   
  ARG=['y             ';
       'rango_retardos';
       'IA            ';];
  
  for i=nargin+1:size(ARG,1)
    eval(sprintf('%s=[];',ARG(i,:)))
  end

  %% VALORES POR DEFECTO 
  if isempty(rango_retardos); rango_retardos=1:(T/3); end 
  if isempty(IA), IA=0; end

  rango_retardos=sort(rango_retardos);
  %% C�LCULO DEL MODELO AR 
  if size(rango_retardos)==[1 1] % si se fija el n�mero de retardos 
    n_retardos=rango_retardos;
    [AR,s2_e]=modar(y,rango_retardos); % Modelo Autorregresivo        
  else				% si se da un rango de retardos 
    j=0; 
    for p=rango_retardos
      j=j+1; 
      [AR,s2_e]=modar(y,p);	% estima el Modelo Autorregresivo 
      BIC(j,1)=log(s2_e)+p*log(T-p)/(T-p); % Crit. de Inf. Bayesiano(BJ,pp.201)
    end 
    [aic,posicion]=min(BIC);
    
    n_retardos=rango_retardos(posicion); % n� retardos �ptimo
    [AR,s2_e]=modar(y,n_retardos); % estima el Modelo Autorregres. 
  end     

  %% INFORMACI�N ADICIONAL 
  if IA>1 
    fprintf('\n\tAR(%1.f)\n\n',n_retardos); 
    if length(rango_retardos)>1 % si se da un rango de retardos 
      plot(rango_retardos,BIC);
      title(sprintf('Criterio de Informaci�n Bayesiano. M�nimo en %g retardos',... 
       n_retardos))
      xlabel('N�mero de retardos')
      ylabel('Valor de la funci�n BIC')
    end 
    if n_retardos>floor(T/3) 
      disp('N�mero de retardos mayor que el 33% del tama�o de la muestra') 
    end 
  end   
