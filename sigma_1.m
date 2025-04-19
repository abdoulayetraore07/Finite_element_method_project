function val = sigma_1(x,y,probleme_aborde)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sigma_1 :
% Evaluation de la fonction sigma_1.
%
% SYNOPSIS val = sigma_1(x,y)
%
% INPUT * x,y : les 2 coordonnees du point ou on veut evaluer la fonction.
%
% OUTPUT - val: valeur de la fonction sur ce point.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(probleme_aborde, 'validation')
    val = 1;

elseif strcmp(probleme_aborde, 'pb_stationnaire_dirichlet_homogene_cas_1')
    val = 5;

elseif strcmp(probleme_aborde, 'pb_stationnaire_dirichlet_homogene_cas_2')
    val = 5;

elseif strcmp(probleme_aborde, 'pb_stationnaire_dirichlet_non_homogene_cas_1')
    val = 1;

elseif strcmp(probleme_aborde, 'pb_stationnaire_dirichlet_non_homogene_cas_2')
    val = 1;

elseif strcmp(probleme_aborde, 'pb_stationnaire_fourier_validation')
    val = 1;

elseif strcmp(probleme_aborde, 'pb_stationnaire_fourier_cas_1')
    val = 1;

elseif strcmp(probleme_aborde, 'pb_stationnaire_fourier_cas_2')
    val = 5;

elseif strcmp(probleme_aborde, 'pb_temporel')
    val = 5;
else
    val= 8;
endif




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%2024%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
