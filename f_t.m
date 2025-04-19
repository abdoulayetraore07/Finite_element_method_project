function val = f_t(x,y,t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f_t :
% Evaluation de la fonction second membre.
%
% SYNOPSIS val = f_t(x,y,t)
%
% INPUT * x,y : les 2 coordonnees du point ou on veut evaluer la fonction.
%       *   t : le temps
%
% OUTPUT - val: valeur de la fonction sur ce point.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

val = exp(-5*t)*600 * exp(-((x - 5)^2 / (0.8^2)) - ((y - 4)^2 / (0.8^2))) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%2024%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
