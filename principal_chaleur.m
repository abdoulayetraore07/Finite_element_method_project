% =====================================================
%
% principal_chaleur;
%
% une routine pour la mise en oeuvre des EF P1 Lagrange
% pour
% 1) l'equation de la chaleur suivante stationnaire, avec condition de
% Dirichlet non homogene
%
% | \alpha T - div(\sigma \grad T)= S,   dans \Omega=\Omega_1 U \Omega_2
% |         T = T_\Gamma,   sur le bord
%
% ou S est la source de chaleur, T_\Gamma la temperature exterieure
% \alpha > 0 et
% \sigma = | \sigma_1 dans \Omega_1
%          | \sigma_2 dans \Omega_2
%
% 2) l'equation de la chaleur dependant du temps avec condition de
% Dirichlet non homogene
%
% | dT/dt - div(\sigma \grad T)= S,   dans \Omega=\Omega_1 U \Omega_2 et pour tout t< t_max
% |         T = T_\Gamma,   sur le bord et pour tout t< t_max
% |         T = T_0       dans \Omega et pour t=0
%
% ou S est la source de chaleur, T_\Gamma la temperature exterieure,
% T_0 est la valeur initiale de la temp?rature
% \alpha > 0 et
% \sigma = | \sigma_1 dans \Omega_1
%          | \sigma_2 dans \Omega_2
% =====================================================


% METTRE LA VALEUR A 'oui' POUR LE PROBLEME DONT ON SOUHAITE UNE SIMULATION ( AUCUNE AUTRE ACTION NECESSAIRE )
% Le pas du maillage est par défaut choisi égale à 1.




%Choix probleme avec validation par defaut
validation = 'non';
pb_stationnaire_dirichlet_homogene_cas_1 = 'non';
pb_stationnaire_dirichlet_homogene_cas_2 = 'non';
pb_stationnaire_dirichlet_non_homogene_cas_1 = 'non';
pb_stationnaire_dirichlet_non_homogene_cas_2 = 'non';
pb_stationnaire_fourier_validation = 'non';
pb_stationnaire_fourier_cas_1 = 'non';
pb_stationnaire_fourier_cas_2 = 'non';
pb_temporel = 'oui';




% Sélection du cas et des paramètres associés
if strcmp(validation, 'oui')
    alpha = 1;
    probleme_aborde = 'validation';
    nom_maillage_msh = 'geomRectangle.msh';
    nom_maillage_geo = 'geomRectangle.geo';

elseif strcmp(pb_stationnaire_dirichlet_homogene_cas_1, 'oui')
    alpha = 1;
    probleme_aborde = 'pb_stationnaire_dirichlet_homogene_cas_1';
    nom_maillage_msh = 'domaine.msh';
    nom_maillage_geo = 'domaine.geo' ;

elseif strcmp(pb_stationnaire_dirichlet_homogene_cas_2, 'oui')
    alpha = 1;
    probleme_aborde = 'pb_stationnaire_dirichlet_homogene_cas_2';
    nom_maillage_msh = 'domaine.msh';
    nom_maillage_geo = 'domaine.geo' ;

elseif strcmp(pb_stationnaire_dirichlet_non_homogene_cas_1, 'oui')
    alpha = 1;
    t_Gamma = 40;
    probleme_aborde = 'pb_stationnaire_dirichlet_non_homogene_cas_1';
    nom_maillage_msh = 'domaine.msh';
    nom_maillage_geo = 'domaine.geo' ;

elseif strcmp(pb_stationnaire_dirichlet_non_homogene_cas_2, 'oui')
    alpha = 1.8;
    t_Gamma = 50;
    probleme_aborde = 'pb_stationnaire_dirichlet_non_homogene_cas_2';
    nom_maillage_msh = 'domaine.msh';
    nom_maillage_geo = 'domaine.geo' ;

elseif strcmp(pb_stationnaire_fourier_validation, 'oui')
    alpha = 1;
    lambda=10 ;
    probleme_aborde = 'pb_stationnaire_fourier_validation';
    nom_maillage_msh = 'domaine.msh';
    nom_maillage_geo = 'domaine.geo' ;

elseif strcmp(pb_stationnaire_fourier_cas_1, 'oui')
    alpha = 2;
    lambda= 0;
    probleme_aborde = 'pb_stationnaire_fourier_cas_1';
    nom_maillage_msh = 'geomRectangle.msh';
    nom_maillage_geo = 'geomRectangle.geo';

elseif strcmp(pb_stationnaire_fourier_cas_2, 'oui')
    alpha = 1;
    lambda= 1;
    probleme_aborde = 'pb_stationnaire_fourier_cas_2';
    nom_maillage_msh = 'domaine.msh';
    nom_maillage_geo = 'domaine.geo' ;

elseif strcmp(pb_temporel, 'oui')
    Tps_initial = 0;
    Tps_final = 1;
    lambda= 1;
    delta_t = 0.01;
    alpha = 1/delta_t;
    N_t = (Tps_final-Tps_initial)/delta_t; % le nombre d'iterations necessaires
    t_Gamma = 0;
    probleme_aborde='pb_temporel';
    nom_maillage_msh = 'domaine.msh';
    nom_maillage_geo = 'domaine.geo' ;

else
    error('Aucun cas sélectionné. Veuillez définir un cas à "oui".');
end



%Affichage du probleme et du ficher de maillage associé
disp(['Cas sélectionné : ', probleme_aborde]);
disp(['Fichier de maillage utilisé : ', nom_maillage_geo]);





% Donnees du probleme et affichage maillage
h = 0.1;
system(['gmsh -2 -clmax ' num2str(h) ' -clmin ' num2str(h) nom_maillage_geo ]);
nom_maillage = nom_maillage_msh;
titre_1 = ['Maillage avec h=' num2str(h) ' et fichier=' nom_maillage];
affichemaillage(nom_maillage, titre_1);
titre_2 = ['Solution u avec h=' num2str(h) ' et maillage=' nom_maillage];





% lecture du maillage et affichage
% ---------------------------------
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(nom_maillage);

% ----------------------
% calcul des matrices EF
% ----------------------

% declarations
% ------------
KK = sparse(Nbpt,Nbpt); % matrice de rigidite
MM = sparse(Nbpt,Nbpt); % matrice de rigidite
LL = zeros(Nbpt,1);     % vecteur second membre
FF = zeros(Nbpt,1);     % vecteur second membre cas fourier_stationnaire
S = zeros(Nbpt,1);     % composant vecteur second membre cas fourier_stationnaire
Uc = zeros(Nbpt,1);     % composant vecteur second membre cas fourier_stationnaire
SS = sparse(Nbpt,Nbpt); % matrice de masse surfacique
t_gamma = zeros(Nbpt,1); % vecteur condtions de bord pour les 3 premiers cas


% boucle sur les triangles
% ------------------------

for l=1:Nbtri

  % calcul des matrices elementaires du triangle l

   [Kel]=matK_elem(Coorneu(Numtri(l,1),:),Coorneu(Numtri(l,2),:),Coorneu(Numtri(l,3),:),Reftri(l),probleme_aborde);

   [Mel]=matM_elem(Coorneu(Numtri(l,1),:),Coorneu(Numtri(l,2),:),Coorneu(Numtri(l,3),:));

    % On fait l'assemblage de la matrice globale

    for i=1:3
       I = Numtri(l, i);
            for j = 1:3
                J = Numtri(l, j);
                MM(I, J) = MM(I, J) + Mel(i, j);
                KK(I, J) = KK(I, J) + Kel(i, j);
            end
    end % for i

end ; % for l


% boucle sur les aretes
% ------------------------

if (strcmp(probleme_aborde,'pb_stationnaire_fourier_validation') || strcmp(probleme_aborde,'pb_stationnaire_fourier_cas_1') || strcmp(probleme_aborde,'pb_stationnaire_fourier_cas_1'))
  for l=1:Nbaretes

      ref=Refaretes(l);
    % calcul des matrices de masse surfacique elementaires de l'arete l
      [Sel] = mat_elem_surface(Coorneu(Numaretes(l,1),:),Coorneu(Numaretes(l,2),:), ref);

      % On fait l'assemblage de la matrice globale

      for i=1:2
         I = Numaretes(l, i);
              for j = 1:2
                  J = Numaretes(l, j);
                  SS(I, J) = SS(I, J) + Sel(i, j);
              end % for j
      end % for i

  end ; % for l
endif



% Matrice EF
% -------------------------
if (strcmp(probleme_aborde,'pb_stationnaire_fourier_validation') || strcmp(probleme_aborde,'pb_stationnaire_fourier_cas_1') || strcmp(probleme_aborde,'pb_stationnaire_fourier_cas_2'))
  AA = alpha*MM+KK+lambda*SS ;
elseif strcmp(probleme_aborde,'pb_temporel')
  AA = (1/delta_t)*MM + KK ;
else
  AA = alpha*MM+KK ;
endif





% =====================================================
% Pour les problemes stationnaires
% ------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcul du second membre L ou FF


if (strcmp(probleme_aborde, 'validation') || strcmp(probleme_aborde, 'pb_stationnaire_dirichlet_homogene_cas_1') || strcmp(probleme_aborde, 'pb_stationnaire_dirichlet_homogene_cas_2'))
    reftri = 1;
    refneu = 1;
    for i = 1:Nbpt
        [val_1, val_2] = f(Coorneu(i, 1), Coorneu(i, 2), probleme_aborde, reftri, refneu);
        LL(i) = val_1;
    end
    LL = MM * LL;

elseif (strcmp(probleme_aborde, 'pb_stationnaire_dirichlet_non_homogene_cas_1'))
    for i = 1:Nbpt
        refneu = 1;
        j = 0; % Identifier dans quelle zone omega 1 ou 2 se trouve un sommet
        for l = 1:Nbtri
            if ((Numtri(l, 1) == i) || (Numtri(l, 2) == i) || (Numtri(l, 3) == i))
                j = l;
                break;
            end
        end
        [val_1, val_2] = f(Coorneu(i, 1), Coorneu(i, 2), probleme_aborde, Reftri(j), refneu);
        LL(i) = val_1;
    end
    LL = MM * LL;

elseif (strcmp(probleme_aborde, 'pb_stationnaire_dirichlet_non_homogene_cas_2'))
    for i = 1:Nbpt
        reftri = 1;
        refneu = 1;
        [val_1, val_2] = f(Coorneu(i, 1), Coorneu(i, 2), probleme_aborde, reftri, refneu);
        LL(i) = val_1;
    end
    LL = MM * LL;

else
    for i = 1:Nbpt
        reftri = 1;
        [val_1, val_2] = f(Coorneu(i, 1), Coorneu(i, 2), probleme_aborde, reftri, Refneu(i));
        S(i, 1) = val_1;
        Uc(i, 1) = val_2;
    end
    FF = MM * S + lambda * SS * Uc;
end










%TECHNIQUE DE PSEUDO_ELIMINATION via fonction elimine.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inversion
% tilde_AA ET tilde_LL SONT LA MATRICE EF ET LE VECTEUR SECOND MEMBRE
% APRES PSEUDO_ELIMINATION

if (strcmp(probleme_aborde, 'validation') || strcmp(probleme_aborde, 'pb_stationnaire_dirichlet_homogene_cas_1') || strcmp(probleme_aborde, 'pb_stationnaire_dirichlet_homogene_cas_2'))
    [tilde_AA, tilde_LL] = elimine(AA, LL, Refneu, Coorneu, probleme_aborde);
    UU = tilde_AA \ tilde_LL;
    for i = 1:Nbpt %%% Ajouter les valeurs aux bords de t_gamma
        %if (Refneu(i) == 1)
        t_gamma(i, 1) = T_Gamma(Coorneu(i, 1), Coorneu(i, 2), probleme_aborde);
        %end
    end
    TT = UU + t_gamma;

elseif (strcmp(probleme_aborde, 'pb_stationnaire_dirichlet_non_homogene_cas_1') || strcmp(probleme_aborde, 'pb_stationnaire_dirichlet_non_homogene_cas_2'))
    [tilde_AA, tilde_LL] = elimine(AA, LL, Refneu, Coorneu, probleme_aborde);
    UU = tilde_AA \ tilde_LL;
    TT = UU;

elseif (strcmp(probleme_aborde, 'pb_stationnaire_fourier_validation') || strcmp(probleme_aborde, 'pb_stationnaire_fourier_cas_1') || strcmp(probleme_aborde, 'pb_stationnaire_fourier_cas_2'))
    UU = AA \ FF;
    TT = UU;

else
    UU = AA \ FF;
    TT = UU;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Calcul des valeurs max et min de TT ( temperature )

[valeur_max, indice_max] = max(TT);
[valeur_min, indice_min] = min(TT);

% Afficher les résultats
disp(['La valeur maximale de T_h est : ', num2str(valeur_max)]);
disp(['La valeur minimale de T_h est : ', num2str(valeur_min)]);







%CALCULS ERREURS EN NORME L2 ET SEMI-NORME H1 pour le cas 'validation'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (strcmp('validation', probleme_aborde))
    TT_exact = 2 * sin(2 * pi * Coorneu(:, 1)) .* sin(3 * pi * Coorneu(:, 2));
    % Calcul de l'erreur L2
    BB = TT_exact - TT;
    err_L2 = sqrt(transpose(BB) * MM * BB);
    % Calcul de l'erreur H1
    err_H1 = sqrt(transpose(BB) * KK * BB);

    % Affichage des erreurs
    disp(['Erreur L2 : ', num2str(err_L2)]);
    disp(['Erreur H1 : ', num2str(err_H1)]);

    % Valeurs de h
    h = [1, 0.4, 0.2, 0.1, 0.05];

    % Erreurs en norme L2 et en semi-norme H1 correspondantes
    L2_errors = [4.8444, 2.1316, 0.40147, 0.14216, 0.038383];
    H1_errors = [23.3665, 7.003, 4.8055, 1.7224, 0.516];

    % Calcul des erreurs relatives
    L2_rel_errors = L2_errors / norm(L2_errors);
    H1_rel_errors = H1_errors / norm(H1_errors);

    % Tracé graphe d'erreurs en norme L2 et en semi-norme H1
    figure;
    loglog(1 ./ h, L2_rel_errors, 'o-', 'LineWidth', 1.5); % Points pour la norme L2
    hold on;

    % Ajustement polynomial pour L2
    p_L2 = polyfit(log(1 ./ h), log(L2_rel_errors), 1);
    loglog(1 ./ h, exp(polyval(p_L2, log(1 ./ h))), 'b--', 'LineWidth', 1.5); % Droite de régression pour L2

    % Points pour la semi-norme H1
    loglog(1 ./ h, H1_rel_errors, 's-', 'LineWidth', 1.5);

    % Ajustement polynomial pour H1
    p_H1 = polyfit(log(1 ./ h), log(H1_rel_errors), 1);
    loglog(1 ./ h, exp(polyval(p_H1, log(1 ./ h))), 'r--', 'LineWidth', 1.5); % Droite de régression pour H1

    % Légende et mise en forme
    legend('Erreur L^2', sprintf('L2 polyfit, pente = %.2f', p_L2(1)), 'Erreur semi-norme H^1', sprintf('H1 polyfit, pente = %.2f', p_H1(1)));
    xlabel('log(1/h)');
    ylabel('log(erreur relative)');
    title('Convergence des erreurs L^2 et semi-norme H^1');
    grid on;
    hold off;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







% Pour le probleme temporel
% ---------------------------------
if (strcmp('pb_temporel',probleme_aborde))

    % on initialise la condition initiale
    % -----------------------------------
    T_initial = condition_initiale(Coorneu(:,1),Coorneu(:,2));

	% solution a t=0
	% --------------
    S = zeros(Nbpt,1);
    TT = 310*ones(Nbpt,1);
    UU = 310*ones(Nbpt,1);

    % visualisation
    % -------------
    figure;
    hold on;
    mot_affichage = strcat(probleme_aborde, ' - %s');
    affiche(TT, Numtri, Coorneu, sprintf(mot_affichage, titre_2));
    %affiche(TT, Numtri, Coorneu, ['Temps = ', num2str(0)]);
    axis([min(Coorneu(:,1)),max(Coorneu(:,1)),min(Coorneu(:,2)),max(Coorneu(:,2)),...
        290,330,290,300]);
    hold off;

	% Boucle sur les pas de temps
	% ---------------------------
    for k = 1:N_t

        LL_k = zeros(Nbpt,1);
        S = zeros(Nbpt,1);

        for i = 1:Nbpt
            x=Coorneu(i,1);
            y=Coorneu(i,2);
            S(i,1) =  f_t(x,y,k);
        end

    % Calcul du second membre F a l instant k*delta t
    % -----------------------------------------------
		% A COMPLETER EN UTILISANT LA ROUTINE f_t.m et le terme precedent (donne par UU)

		LL_k =    (1/delta_t)*MM*UU + MM*S;

		% inversion
		% ----------
		% tilde_AA ET tilde_LL_k SONT LA MATRICE EF ET LE VECTEUR SECOND MEMBRE
		% APRES PSEUDO_ELIMINATION
		% ECRIRE LA ROUTINE elimine.m ET INSERER L APPEL A CETTE ROUTINE
		% A UN ENDROIT APPROPRIE
    [tilde_AA, tilde_LL_k] = elimine(AA, LL_k, Refneu, Coorneu, probleme_aborde);
    UU = tilde_AA\tilde_LL_k;
    TT = UU ;

        % visualisation
		%& -------------
        pause(0.05)
        mot_affichage = strcat(probleme_aborde, ' - %s');
        %affiche(TT, Numtri, Coorneu, ['Temps = ', num2str(k*delta_t)]);
        affiche(TT, Numtri, Coorneu, sprintf(mot_affichage, titre_2));
        axis([min(Coorneu(:,1)),max(Coorneu(:,1)),min(Coorneu(:,2)),max(Coorneu(:,2)),...
            290,330,290,320]);
    end

end










% visualisation

% Création de la chaîne de caractères pour l'affichage


%Enfin, l'affichage
if (strcmp('pb_temporel',probleme_aborde))==0
  affiche(TT, Numtri, Coorneu, sprintf(mot_affichage, titre_2));
end

