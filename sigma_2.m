function val = sigma_2(x,y,probleme_aborde)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sigma_2 :
% Evaluation de la fonction sigma_2.
%
% SYNOPSIS val = sigma_2(x,y)
%
% INPUT * x,y : les 2 coordonnees du point ou on veut evaluer la fonction.
%
% OUTPUT - val: valeur de la fonction sur ce point.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(probleme_aborde, 'validation')
    val = 1;

elseif strcmp(probleme_aborde, 'pb_stationnaire_dirichlet_homogene_cas_1')
    val = sqrt(3)/2;

elseif strcmp(probleme_aborde, 'pb_stationnaire_dirichlet_homogene_cas_2')
    val = 1/4*( (2+sin(16*pi*x))*(2+sin(16*pi*y))  );

elseif strcmp(probleme_aborde, 'pb_stationnaire_dirichlet_non_homogene_cas_1')
    val = sqrt(3)/2;

elseif strcmp(probleme_aborde, 'pb_stationnaire_dirichlet_non_homogene_cas_2')
    val = sqrt(3)/2;

elseif strcmp(probleme_aborde, 'pb_stationnaire_fourier_validation')
    val = 0.5;

elseif strcmp(probleme_aborde, 'pb_stationnaire_fourier_cas_1')
    val = 1;

elseif strcmp(probleme_aborde, 'pb_stationnaire_fourier_cas_2')
    val =  1/4*( (2+sin(16*pi*x))*(2+sin(16*pi*y)) );

elseif strcmp(probleme_aborde, 'pb_temporel')
    val = 1/4*( (2+sin(16*pi*x))*(2+sin(16*pi*y))  );
else
    val = 0; % Cas par défaut
endif


% A CHANGER POUR LA VALIDATION


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%2024%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
