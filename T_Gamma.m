function val = T_Gamma(x,y,probleme_aborde)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T_Gamma :
% renvoie la valeur de T_Gamma,h a un sommet dans la routine elimine
%
% SYNOPSIS val = T_Gamma(S)
%
% INPUT * S : numero d'un sommet
%
% OUTPUT - val : valeur de T_Gamma,h au sommet numero S.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%val= sin(pi*x)*cos(pi*y) ;
%if y==6
%  val=;
%else
%  val=0;
%endif

if strcmp(probleme_aborde, 'validation')
    val = 0;

elseif strcmp(probleme_aborde, 'pb_stationnaire_dirichlet_homogene_cas_1')
    val = 285;

elseif strcmp(probleme_aborde, 'pb_stationnaire_dirichlet_homogene_cas_2')
    val = 285;

elseif strcmp(probleme_aborde, 'pb_stationnaire_dirichlet_non_homogene_cas_1')
    val = sin(pi*x)*cos(pi*y) ;

elseif strcmp(probleme_aborde, 'pb_stationnaire_dirichlet_non_homogene_cas_2')
    if y==6
      val=310;
    else
      val=260;
    endif
else
    val= 285;
endif


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
