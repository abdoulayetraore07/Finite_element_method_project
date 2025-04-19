function [Sel] = mat_elem_surface(S1, S2, ref)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mat_elem :
% calcul les matrices de masse surfaciques pour les frontieres
% Gamma1, Gamma2, Gamma3 et Gamma4
%
% SYNOPSIS [Sel] = mat_elem_surface(S1, S2, ref)
%
% INPUT * S1, S2 : les 2 coordonnees des 2 sommets de l'arete
%                      (vecteurs reels 1x2)
%         ref    : Reference de l'arete.
%
% OUTPUT - Sel matrice de masse surfacique elementaire pour le bord de
% reference ref
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% preliminaires, pour faciliter la lecture:
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);

% On calcul la longueur de l'arete.
Long = sqrt((x2-x1)^2+(y2-y1)^2);

% Declaration et ramplisage des matrice elementaires.
%---------------------------------------------------
Sel = zeros(2,2);
if ref==1
  for i=1:2
	  for j=1:2
		  if i==j
        Sel(i,j) = (1/12)*abs(Long) ;
      else
        Sel(i,j) = (1/24)*abs(Long) ;
      endif
	  end; % j
  end; % i
endif

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%2024%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
