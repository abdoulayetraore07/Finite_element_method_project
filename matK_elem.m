function [Kel] = matK_elem(S1, S2, S3, ref,probleme_aborde)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mat_elem :
% calcul la matrice de raideur élémentaire en P1 Lagrange
%
% SYNOPSIS [Kel] = mat_elem(S1, S2, S3)
%
% INPUT * S1, S2, S3 : les coordonnées des 3 sommets du triangle
%                      (vecteurs réels 1x2)
%
% OUTPUT - Kel matrice de raideur élémentaire (matrice 3x3)
%
% NOTE (1) Utilisation d'une quadrature à 3 points d'ordre 2
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% préliminaires, pour faciliter la lecture :
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

% les 3 normales à l'arête opposée (de la longueur de l'arête)
norm = zeros(3, 2);
norm(1, :) = [y2-y3, x3-x2];
norm(2, :) = [y3-y1, x1-x3];
norm(3, :) = [y1-y2, x2-x1];

% D est, au signe près, deux fois l'aire du triangle
D = ((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1));
if (abs(D) <= eps)
  error('L''aire d''un triangle est nulle!!!');
end;

% définition des matrices Bl et Sl
Bl = zeros(2,2); Sl = zeros(2,1);
Bl(1,1) = (x2-x1); Bl(1,2) = x3-x1; Bl(2,1) = y2-y1; Bl(2,2) = y3-y1;
Sl(1,1) = x1; Sl(2,1) = y1;

% définition de wo, points de quadrature et gradients des w_i avec grad w_i en colonne i ;
wo = 1/6; so = 1/6; s1 = 2/3;
grad_w = zeros(2,3);
grad_w(1,1) = -1; grad_w(2,1) = -1;
grad_w(1,2) = 1; grad_w(2,2) = 0;
grad_w(1,3) = 0; grad_w(2,3) = 1;

% définition de la fonction G_i,j
function G = G(i,j,x,y)
  M = zeros(2,1); M(1,1) = x; M(2,1) = y;
  Fl = Bl * M + Sl;
  if ref == 1
    sigma_Fl = sigma_1(Fl(1,1), Fl(2,1),probleme_aborde);
  else
    sigma_Fl = sigma_2(Fl(1,1), Fl(2,1),probleme_aborde);
  end
  v1 = transpose(Bl) \ grad_w(:,i); % Résout Bl' * v1 = grad_w(:,i)
  v2 = transpose(Bl) \ grad_w(:,j); % Résout Bl' * v2 = grad_w(:,j)
  G = sigma_Fl * dot(v1, v2);

end

% calcul de la matrice de raideur
% -------------------------------
Kel = zeros(3,3);
for i = 1:3
  for j = 1:3
    Kel(i,j) = abs(D) * wo * ( G(i,j,so,so) + G(i,j,s1,so) + G(i,j,so,s1) )
  end
end;

end

