function [tilde_AA, tilde_LL] = elimine(AA, LL, Refneu, Coorneu, probleme_aborde)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% elimine :
% Realise la pseudo-élimination des nœuds du bord en utilisant le vecteur Refneu.
%
% SYNOPSIS [tilde_AA, tilde_LL] = elimine(AA, LL, Refneu, Coorneu, probleme_aborde)
%
% INPUT * AA, LL, Refneu, Coorneu, probleme_aborde :
%                          la matrice AA assemblée avant élimination,
%                          La matrice LL second membre
%                          le tableau Refneu reference des sommets (vecteur entier Nbpt x 1)
%                          Coorneu
%                          probleme_aborde: pour spécifier le probleme abordé et eviter de lourds petits ajustements
%
% OUTPUT - tilde_AA, tilde_LL : matrice assemblée et second membre
%                               après éléimination des noeuds
%                               du bord
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(probleme_aborde, 'validation')
  for i=1:length(LL)
    if  Refneu(i)==1
      AA(i,:) = 0 ;  AA(:,i) = 0  ;
      LL(i) = 0 ;
    endif
  end % for i

else
  for i=1:length(LL)
    if  Refneu(i)==1
      AA(i,:) = 0 ;
      x=Coorneu(i,1); y=Coorneu(i,2);
      LL(i) = T_Gamma(x,y,probleme_aborde) ;
    endif
  end % for i
end



tilde_AA = AA ; tilde_LL = LL ;


