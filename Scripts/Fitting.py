import numpy as np
import scipy.linalg

#=======================================================================================================================

# Renvoit les coordonnées des points de la surface moyenne linéaire ou quadratique d'un nuage de point en 3D
## x, y, z : coordonnées des points du nuage selon les axes x, y et z respectivement
## order : 1 pour linéaire et 2 pour quadratique
### X, Y, Z : coordonnées des points du plan moyen du nuage selon les axes x, y et z respectivement

def plane_fitting(x, y, z, order):
    data = np.c_[x, y, z]
    X, Z = np.meshgrid(np.arange(min(x), max(x), 0.5), np.arange(min(z), max(z), 0.5))
    XX = X.flatten()
    ZZ = Z.flatten()
    if order == 1:
        A = np.c_[data[:, 0], data[:, 2], np.ones(data.shape[0])]
        C, _, _, _ = scipy.linalg.lstsq(A, data[:, 1])
        Y = C[0] * X + C[1] * Z + C[2]
    elif order == 2:
        A = np.c_[np.ones(data.shape[0]), data[:, [0, 2]], np.prod(data[:, [0, 2]], axis=1), data[:, [0, 2]] ** 2]
        C, _, _, _ = scipy.linalg.lstsq(A, data[:, 1])
        Y = np.dot(np.c_[np.ones(XX.shape), XX, ZZ, XX * ZZ, XX ** 2, ZZ ** 2], C).reshape(X.shape)
    return X, Y, Z

#=======================================================================================================================

# Renvoit les coordonnées des points de la courbe moyenne linéaire ou quadratique d'un nuage de point en 2D
## x, y : coordonnées des points du nuage selon les axes x et y respectivement
## order : 1 pour linéaire et 2 pour quadratique
### courbe : liste des coefficients caractéristiques de la courbe par degré croissant

def line_fitting(x, y, order):
    courbe = np.poly1d(np.polyfit(x, y, order))
    return courbe

