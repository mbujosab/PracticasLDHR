%
% SIDINIT  - Initializes the default options stored in the global vector
%            SIDOPTION
%
% This function should be executed at least once at the beginning of a session
%
% 23/12/96
% (C@@)

global SIDOPTION SIDERROR SIDWARN SIDDISP
if isempty(SIDOPTION)
   SIDOPTION = zeros(1,10);
end
sidopt;

% Inicializaci�n de los mensajes de error
SIDERROR = ['1. Ejecute SIDINIT antes de usar la librer�a                         ';
            '2. N�mero de argumentos incorrecto                                   ';
            '3. SIDOPT: Opci�n no reconocida: %s                                  ';
            '4. SIDOPT: Valor no reconocido: %s                                   ';
            '5. n debe ser menor o igual a i y a j                                '];
            
% Inicializaci�n de los mensajes de aviso
SIDWARN = [];

% Inicializaci�n de mensajes de informaci�n
SIDDISP = [];

