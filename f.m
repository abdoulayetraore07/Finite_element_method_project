function [val_1,val_2] = f(x,y, probleme_aborde, Reftri,Refneu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f :
% Evaluation de la fonction second membre.
%
% SYNOPSIS val = f(x,y)
%
% INPUT * x,y : les 2 coordonnees du point ou on veut evaluer la fonction.
%
% OUTPUT - val: valeur de la fonction sur ce point.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(probleme_aborde, 'validation')
    val_1 = (2 + 26 * pi^2) * sin(2 * pi * x) * sin(3 * pi * y);
    val_2=0;

elseif strcmp(probleme_aborde, 'pb_stationnaire_dirichlet_homogene_cas_1')
     val_1 = 600 * exp(-((x - 5)^2 / (0.8^2)) - ((y - 4)^2 / (0.8^2))) - 285;
     val_2=0 ;

elseif strcmp(probleme_aborde, 'pb_stationnaire_dirichlet_homogene_cas_2')
     val_1 = 600 * exp(-((x - 5)^2 / (0.8^2)) - ((y - 4)^2 / (0.8^2))) - 285;
     val_2=0 ;

elseif strcmp(probleme_aborde, 'pb_stationnaire_dirichlet_non_homogene_cas_1')
  if Reftri==1
      val=(1 + 2 * pi^2) * sin(pi * x) * cos(pi * y);
      val=0;
  else
      val=(1 + sqrt(3) * pi^2) * sin(pi * x) * cos(pi * y);
      val=0;
  endif

elseif strcmp(probleme_aborde, 'pb_stationnaire_dirichlet_non_homogene_cas_2')
    val_1 = 0;
    val_2= 0;

elseif strcmp(probleme_aborde, 'pb_stationnaire_fourier_validation')
    val_1= 0;
    val_2= 0;
    if (Refneu==1) &&  ( (x==0) && (y>=0) && (y<=6)  )
      val_2= 20 ;
    endif

elseif strcmp(probleme_aborde, 'pb_stationnaire_fourier_cas_1')
    val_1 = sin(pi*x/2)*sin(pi*y/2);
    val_2 = 0;

elseif strcmp(probleme_aborde, 'pb_stationnaire_fourier_cas_2')
    val_1 = 600 * exp(-((x - 5)^2 / (0.8^2)) - ((y - 4)^2 / (0.8^2))) ;
    val_2 = 0 ;
    if Refneu==1
      val_2 = 285 ;
    endif

elseif strcmp(probleme_aborde, 'pb_temporel')
    val_1 = 9;
    val_2=0;

else
    val_1= 8;
    val_2=0;

endif
% A CHANGER POUR LA VALIDATION


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%2024%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
