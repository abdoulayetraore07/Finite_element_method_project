Mesh.MshFileVersion = 2.2;
// definition du pas du maillage
h = 0.1;
// définition des points (en 3D, raison pour laquelle il y a un 0 en z)
Point(1) = {0, 0, 0, h};
Point(2) = {2, 0, 0, h};
Point(3) = {2, 2, 0, h};
Point(4) = {0, 2, 0, h};
Point(5) = {0.2, 0.2, 0, h};
Point(6) = {1.8, 0.2, 0, h};
Point(7) = {1.8, 1.2, 0, h};
Point(8) = {0.2, 1.2, 0, h};
// définition des segments qui relient les points
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
// définition des contours fermés
Line Loop(1) = {1,2,3,4,5,6,7,8};
Line Loop(2) = {5,6,7,8};
// définition des surfaces à partir contours fermés
Plane Surface(1) = {1};
Plane Surface(2) = {2};
// définition des éléments physiques : pour ces éléments, nous pourrons récupérer
//									   les références 
Physical Point(1) = {1,2,3,4};
Physical Line(1) = {1,2,3,4};
Physical Line(2) = {5,6,7,8};
Physical Surface(1) = {1};
Physical Surface(2) = {2};
