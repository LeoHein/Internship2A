## Ce script présente la facon de créer un catalogue puis un échantillon varié des différents affichages possibles.
# La première partie présente le paramétrage du catalogue
# La seconde partie montre un exemple de catalogue pour L1630S et divers affichages classiques
# La troisième partie montre un exemple de catalogue pour Orion complet et divers affichages classiques
# Le script est à executer en l'état, fermer la figure pour faire apparaître la suivante


# IMPORTS
from Localization import *
from Catalog import *
from MapsPlot import *

#=======================================================================================================================
"""
## PARAMETRAGE ET CREATION DU CATALOGUE
# CHOIX DES PARAMETRES DES DIFFERENTS CATALOGUES PARMI LES LISTES CI-DESSOUS
'angDist, _RAJ2000, _DEJ2000, RAJ2000, DEJ2000, Catalog, OriClass, 3.6mag, e_3.6mag, 4.5mag, e_4.5mag, 5.8mag, e_5.8mag, 8.0mag, e_8.0mag, 24mag, e_24mag, Target, Pred, P(CI), P(CII), P(Other)'
'RAJ2000,  DEJ2000,  2MASS,  Gaia,  more,  Nep,  Teff,  e_Teff,  logg,  e_logg,  l_RVmean,  RVmean,  e_RVmean,  Fbol,  e_Fbol,  Av,  e_Av,  Theta,  e_Theta,  Flag,  Group,  SimbadName'
'angDist, _RAJ2000, _DEJ2000, recno, ID, Nvis, Loc, RAJ2000, DEJ2000, SNR, SaFlag, HRV, e_HRV, s_HRV, errHRV, f_Teff, Teff, e_Teff, f_logg, logg, e_logg, f_Vmicro, Vmicro, Vmacro, f_Vsini, Vsini, f_[M/H], [M/H], e_[M/H], f_[a/M], [a/M], e_[a/M], Chi2, TClass, AFlag, f_[C/Fe], [C/Fe], e_[C/Fe], f_[CI/Fe], [CI/Fe], e_[CI/Fe], f_[N/Fe], [N/Fe], e_[N/Fe], f_[O/Fe], [O/Fe], e_[O/Fe], f_[Na/Fe], [Na/Fe], e_[Na/Fe], f_[Mg/Fe], [Mg/Fe], e_[Mg/Fe], f_[Al/Fe], [Al/Fe], e_[Al/Fe], f_[Si/Fe], [Si/Fe], e_[Si/Fe], f_[P/Fe], [P/Fe], e_[P/Fe], f_[S/Fe], [S/Fe], e_[S/Fe], f_[K/Fe], [K/Fe], e_[K/Fe], f_[Ca/Fe], [Ca/Fe], e_[Ca/Fe], f_[Ti/Fe], [Ti/Fe], e_[Ti/Fe], f_[TiII/Fe], [TiII/Fe], e_[TiII/Fe], f_[V/Fe], [V/Fe], e_[V/Fe], f_[Cr/Fe], [Cr/Fe], e_[Cr/Fe], f_[Mn/Fe], [Mn/Fe], e_[Mn/Fe], f_[Fe/H], [Fe/H], e_[Fe/H], f_[Co/Fe], [Co/Fe], e_[Co/Fe], f_[Ni/Fe], [Ni/Fe], e_[Ni/Fe], f_[Cu/Fe], [Cu/Fe], e_[Cu/Fe], f_[Ce/Fe], [Ce/Fe], e_[Ce/Fe], Gaia, Sp'
'ra_epoch2000, dec_epoch2000, errHalfMaj, errHalfMin, errPosAng, source_id, ra, ra_error, dec, dec_error, parallax, parallax_error, parallax_over_error, pm, pmra, pmra_error, pmdec, pmdec_error, astrometric_n_good_obs_al, astrometric_gof_al, astrometric_chi2_al, astrometric_excess_noise, astrometric_excess_noise_sig, astrometric_params_solved, pseudocolour, pseudocolour_error, visibility_periods_used, ruwe, duplicated_source, phot_g_mean_flux, phot_g_mean_flux_error, phot_g_mean_mag, phot_bp_mean_flux, phot_bp_mean_flux_error, phot_bp_mean_mag, phot_rp_mean_flux, phot_rp_mean_mag, phot_bp_rp_excess_factor, bp_rp, dr2_radial_velocity, dr2_radial_velocity_error, dr2_rv_nb_transits, dr2_rv_template_teff, dr2_rv_template_logg, panstarrs1, sdssdr13, skymapper2, urat1, phot_g_mean_mag_error, phot_bp_mean_mag_error, phot_rp_mean_mag_error, phot_g_mean_mag_corrected, phot_g_mean_mag_error_corrected, phot_g_mean_flux_corrected, phot_bp_rp_excess_factor_corrected, ra_epoch2000_error, dec_epoch2000_error, ra_dec_epoch2000_corr'

param_cornu = ['Pred', 'P(CI)', 'P(CII)']  # Paramètres sélectionnés dans le catalogue de Cornu
param_kounkel = ['2MASS', 'Gaia', 'RVmean', 'e_RVmean', 'Av', 'e_Av']  # Paramètres sélectionnés dans le catalogue de Kounkel
param_johnsson = ['ID', 'HRV', 'errHRV']  # Paramètres sélectionnés dans le catalogue de Jonsson
param_gaia = ['source_id', 'parallax', 'parallax_error', 'pmra', 'pmra_error', 'pmdec', 'pmdec_error']  # Paramètres sélectionnés dans le catalogue de Gaia

# CHOIX DU SEUIL DE PROBABILITE DE CORNU ACCEPTE
ProbThresh = 0  # Pour utiliser la classification de Cornu
ProbThresh = 0.95  # Exemple de seuil

# CHOIX DU CATALOGUE APOGEE
Apogee, Velocity, VelErr = 'Kounkel', 'RVmean', 'e_RVmean'  # Pour utiliser le catalogue de Kounkel et ses vitesses radiales
Apogee, Velocity, VelErr = 'Johnsson', 'HRV', 'errHRV'  # Pour utiliser le catalogue de Jonsson et ses vitesses radiales

# SELECTION DES TYPES D'ETOILES A SUPPRIMER
Remove = 'NomencChange'  # Pour supprimer les étoiles ayant changé de nomenclature entre Gaia DR2 et DR3
Remove = 'double'  # Pour supprimer les étoiles devenues objets doubles entre Gaia DR2 et Gaia DR3
Remove = 'Both'  # Pour supprimer les deux catégories ci-dessus
Remove = 'None'  # Pour ne rien supprimer

# CHOIX DE LA LOCALISATION
Loc = [86,  88,  -1,  1]  # Objets dans L1630N seulement
Loc = [85,  86,  -3,  -0.8]  # Objets dans L1630S seulement
Loc = [83,  86.5,  -11,  -4]  # Objets dans Orion A seulement
Loc = [0,  360,  -90,  90]  # Ensemble du catalogue de Cornu

# CREATION DU CATALOGUE ET EXPORT
Tab = catalog(Param_cornu=param_cornu, Param_Gaia=param_gaia, Param_Johnsson=param_johnsson, Param_Kounkel=param_kounkel, ProbThresh=ProbThresh, Loc=Loc, RemoveType=Remove, Apo=Apogee)
Tab.to_csv('res_cornuXapoXgaia.txt',  header=True,  decimal='.',  index=False,  sep='\t',  float_format='%.7f')

#=======================================================================================================================

## FILTRAGE ET AFFICHAGES DES DONNEES DU CATALOGUE
# PARAMETRAGE DES AFFICHAGES
type = 'Relative'  # Pour des vitesses relatives au mouvement moyen du nuage
type = 0  # Pour des vitesses brutes
abs = 1  # Pour un affichage en valeurs absolues des mouvements propres
abs = 0  # Pour un affichage normal des mouvements propres

# AFFICHAGE DES DONNEES BRUTES
distance_graph(df=Tab)  # Distribution croissante des distances
velocity_graph(df=Tab, velocity=Velocity, velerr=VelErr, type=0)  # Distribution croissante des vitesses radiales
pm_graph(df=Tab,  abs=0)  # Distribution croissante des mouvements propres
velocity_graphs(df=Tab, velocity=Velocity, velerr=VelErr, type=0)  # Distribution croissante des vitesses

# EXEMPLES DE FILTRAGES
# Remarque : les bornes et erreurs sont à changer directement dans l'appel des fonctions
Tab = position_filter(df=Tab, MinDist=300, MaxDist=500, error=50)  # Filtrage des distances
Tab = velocity_filter(df=Tab, velocity=Velocity, velerr=VelErr, vmin=1, vmax=16, error=0.5)  # Filtrage des vitesses radiales
Tab = pm_filter(df=Tab,  pmramin=-10,  pmramax=10,  pmdecmin=-10,  pmdecmax=10,  type='relative',  error=0.5)  # Filtrage des mouvements propres

# AFFICHAGE DES DONNEES FILTREES
distance_graph(df=Tab)  # Distribution croissante des distances
velocity_graph(df=Tab, velocity=Velocity, velerr=VelErr, type=0)  # Distribution croissante des vitesses radiales
pm_graph(df=Tab,  abs=0)  # Distribution croissante des mouvements propres
velocity_graphs(df=Tab, velocity=Velocity, velerr=VelErr, type=0)  # Distribution croissante des vitesses

# EXPORT DU CATALOGUE FILTRE
Tab.to_csv('res_cornuXapoXgaia_filtre.txt',  header=True,  decimal='.',  index=False,  sep='\t',  float_format='%.7f')
"""

#=======================================================================================================================


### EXEMPLE 1 : CATALOGUE POUR L1630S
## CREATION DU CATALOGUE ET FILTRAGE
param_cornu = ['Pred', 'P(CI)', 'P(CII)']
param_kounkel = ['2MASS', 'Gaia', 'RVmean', 'e_RVmean', 'Av', 'e_Av']
param_johnsson = ['ID', 'HRV', 'errHRV']
param_gaia = ['source_id', 'parallax', 'parallax_error', 'pmra', 'pmra_error', 'pmdec', 'pmdec_error']
ProbThresh = 0
Apogee, Velocity, VelErr = 'Kounkel', 'RVmean', 'e_RVmean'
Remove = 'Both'
Loc = [85,  86,  -3,  -0.8]

Tab = catalog(Param_cornu=param_cornu, Param_Gaia=param_gaia, Param_Johnsson=param_johnsson, Param_Kounkel=param_kounkel, ProbThresh=ProbThresh, Loc=Loc, RemoveType=Remove, Apo=Apogee)

Tab = position_filter(df=Tab, MinDist=300, MaxDist=500, error=100)
Tab = velocity_filter(df=Tab, velocity=Velocity, velerr=VelErr, vmin=7, vmax=16, error=1)
Tab = pm_filter(df=Tab,  pmramin=-10,  pmramax=10,  pmdecmin=-10,  pmdecmax=10,  type='relative',  error=0)


## AFFICHAGES 2D DES VITESSES RADIALES
# Remarque : Il faut bien modifier vmin et vmax pour fixer l'echelle de couleur qui dépend de l'échantillon de YSOs
# Remarque : Pour les deux fits, bien stipuler vinf et vsup pour séparer les deux échantillons
# AFFICHAGE NUAGE DE POINTS SIMPLE
velocity_2D(df=Tab, velocity=Velocity, velerr=VelErr, pm=0, fit=0, order=0, vinf=0, vsup=0, vmin=7, vmax=13)
# AFFICHAGE NUAGE DE POINTS AVEC DEUX FITS D'ORDRE 1
velocity_2D(df=Tab, velocity=Velocity, velerr=VelErr, pm=0, fit=2, order=1, vinf=10, vsup=11, vmin=7, vmax=13)
# AFFICHAGE NUAGE DE POINTS AVEC DEUX FITS D'ORDRE 2
velocity_2D(df=Tab, velocity=Velocity, velerr=VelErr, pm=0, fit=2, order=2, vinf=10, vsup=11, vmin=7, vmax=13)
# AFFICHAGE NUAGE DE VECTEURS SIMPLE
velocity_2D(df=Tab, velocity=Velocity, velerr=VelErr, pm=1, fit=0, order=0, vinf=0, vsup=0, vmin=7, vmax=13)
# AFFICHAGE NUAGE DE VECTEURS AVEC DEUX FITS D'ORDRE 1
velocity_2D(df=Tab, velocity=Velocity, velerr=VelErr, pm=1, fit=2, order=1, vinf=10, vsup=11, vmin=7, vmax=13)
# AFFICHAGE NUAGE DE VECTEURS AVEC DEUX FITS D'ORDRE 2
velocity_2D(df=Tab, velocity=Velocity, velerr=VelErr, pm=1, fit=2, order=2, vinf=10, vsup=11, vmin=7, vmax=13)


## AFFICHAGES 3D DES VITESSES
# Remarque : Il faut bien modifier vmin et vmax pour fixer l'echelle de couleur qui dépend de l'échantillon de YSOs
# Remarque : Pour les deux fits, bien stipuler vinf et vsup pour séparer les deux échantillons
# Remarque : Pour le nuage de vecteurs, l'échelle de couleur concerne ici la vitesse totale (cf velocity_graphs)
# Remarque : Pour les deux fits, si le nombre de YSOs n'est pas assez élevé, au moins un ne s'affiche pas
# NUAGE DE POINTS SIMPLE
velocity_3D(df=Tab, velocity=Velocity, velerr=VelErr, pm=0, fit=0, order=0, vinf=0, vsup=0, vmin=7, vmax=13)
# NUAGE DE POINTS AVEC UN FIT D'ORDRE 1
velocity_3D(df=Tab, velocity=Velocity, velerr=VelErr, pm=0, fit=1, order=1, vinf=0, vsup=0, vmin=7, vmax=13)
# NUAGE DE POINTS AVEC UN FIT D'ORDRE 2
velocity_3D(df=Tab, velocity=Velocity, velerr=VelErr, pm=0, fit=1, order=2, vinf=0, vsup=0, vmin=7, vmax=13)
# NUAGE DE POINTS AVEC DEUX FITS D'ORDRE 1
velocity_3D(df=Tab, velocity=Velocity, velerr=VelErr, pm=0, fit=2, order=1, vinf=10.8, vsup=10.8, vmin=7, vmax=13)
# NUAGE DE POINTS AVEC DEUX FITS D'ORDRE 2
velocity_3D(df=Tab, velocity=Velocity, velerr=VelErr, pm=0, fit=2, order=2, vinf=10.8, vsup=10.8, vmin=7, vmax=13)
# NUAGE DE VECTEURS SIMPLE
velocity_3D(df=Tab, velocity=Velocity, velerr=VelErr, pm=1, fit=0, order=0, vinf=0, vsup=0, vmin=1, vmax=7)
# NUAGE DE VECTEURS AVEC UN FIT D'ORDRE 1
velocity_3D(df=Tab, velocity=Velocity, velerr=VelErr, pm=1, fit=1, order=1, vinf=0, vsup=0, vmin=1, vmax=7)
# NUAGE DE VECTEURS AVEC UN FIT D'ORDRE 2
velocity_3D(df=Tab, velocity=Velocity, velerr=VelErr, pm=1, fit=1, order=2, vinf=0, vsup=0, vmin=1, vmax=7)
# NUAGE DE VECTEURS AVEC DEUX FITS D'ORDRE 1
velocity_3D(df=Tab, velocity=Velocity, velerr=VelErr, pm=1, fit=2, order=1, vinf=2, vsup=2, vmin=1, vmax=7)
# NUAGE DE VECTEURS AVEC DEUX FITS D'ORDRE 2
velocity_3D(df=Tab, velocity=Velocity, velerr=VelErr, pm=1, fit=2, order=2, vinf=2, vsup=2, vmin=1, vmax=7)


## EXEMPLES CLASSIQUES DE CARTOGRAPHIES
# NUAGE DE POINTS DE VITESSES RADIALES
mapsplot(df=Tab, plottype='Scatter', maptype='Mean Velocity', region='Main', p='Radial Velocity', velocity=Velocity, zoom=1, axis=[700, 1000, 100, 650], vmin=7, vmax=13)
# NUAGE DE POINTS DE DISTANCES
mapsplot(df=Tab, plottype='Scatter', maptype='Mean Velocity', region='Main', p='Distance', velocity=Velocity, zoom=1, axis=[700, 1000, 100, 650], vmin=300, vmax=450)
# NUAGE DE POINTS D'EXTINCTION
mapsplot(df=Tab, plottype='Scatter', maptype='Mean Velocity', region='Main', p='Extinction', velocity=Velocity, zoom=1, axis=[700, 1000, 100, 650], vmin=10**20, vmax=10**22)
# NUAGE DE POINTS D'EXTINCTION SUR CARTE DE LOMBARDI
mapsplot(df=Tab, plottype='Scatter', maptype='Lombardi', region='Main', p='Extinction', velocity=Velocity, zoom=1, axis=[700, 1000, 100, 650], vmin=10**20, vmax=10**22)
# NUAGE DE POINTS DE DISTANCE SUR CARTE DE LOMBARDI
mapsplot(df=Tab, plottype='Scatter', maptype='Lombardi', region='Main', p='Distance', velocity=Velocity, zoom=1, axis=[700, 1000, 100, 650], vmin=300, vmax=450)
# NUAGE DE POINTS DE VITESSES RADIALES SUR CARTE DE LOMBARDI
mapsplot(df=Tab, plottype='Scatter', maptype='Lombardi', region='Main', p='Radial Velocity', velocity=Velocity, zoom=1, axis=[700, 1000, 100, 650], vmin=7, vmax=13)
# NUAGE DE VECTEURS DE VITESSES RADIALES
mapsplot(df=Tab, plottype='Quiver', maptype='Mean Velocity', region='Main', p='Radial Velocity', velocity=Velocity, zoom=1, axis=[700, 1000, 100, 650], vmin=7, vmax=13)
# NUAGE DE VECTEURS DE DISTANCES
mapsplot(df=Tab, plottype='Quiver', maptype='Mean Velocity', region='Main', p='Distance', velocity=Velocity, zoom=1, axis=[700, 1000, 100, 650], vmin=300, vmax=450)
# NUAGE DE VECTEURS D'EXTINCTION
mapsplot(df=Tab, plottype='Quiver', maptype='Mean Velocity', region='Main', p='Extinction', velocity=Velocity, zoom=1, axis=[700, 1000, 100, 650], vmin=10**20, vmax=10**22)
# NUAGE DE VECTEURS D'EXTINCTION SUR CARTE DE LOMBARDI
mapsplot(df=Tab, plottype='Quiver', maptype='Lombardi', region='Main', p='Extinction', velocity=Velocity, zoom=1, axis=[700, 1000, 100, 650], vmin=10**20, vmax=10**22)
# NUAGE DE VECTEURS DE DISTANCE SUR CARTE DE LOMBARDI
mapsplot(df=Tab, plottype='Quiver', maptype='Lombardi', region='Main', p='Distance', velocity=Velocity, zoom=1, axis=[700, 1000, 100, 650], vmin=300, vmax=450)
# NUAGE DE VECTEURS DE VITESSES RADIALES SUR CARTE DE LOMBARDI
mapsplot(df=Tab, plottype='Quiver', maptype='Lombardi', region='Main', p='Radial Velocity', velocity=Velocity, zoom=1, axis=[700, 1000, 100, 650], vmin=7, vmax=13)


## AFFICHAGE DE CARTES DE DIFFERENCS DE VITESSES
# NUAGE DE POINTS DE DIFFERENCES
comparison(df=Tab, velocity=Velocity, maptype='Mean Velocity', region='Main', plottype='MapDiff', type=0, vmin=-4, vmax=4, mapplot=0, thresh=0)
# NUAGE DE POINTS DE DIFFERENCES EN VALEURS ABSOLUES
comparison(df=Tab, velocity=Velocity, maptype='Mean Velocity', region='Main', plottype='MapDiff', type='abs', vmin=0, vmax=4, mapplot=0, thresh=0)
# NUAGE DE POINTS DE DIFFERENCES AVEC SEUIL LIMITE DE 1.5 AU MINIMUM
comparison(df=Tab, velocity=Velocity, maptype='Mean Velocity', region='Main', plottype='MapDiff', type=0, vmin=-4, vmax=4, mapplot=0, thresh=1.5)
# NUAGE DE POINTS DE DIFFERENCES AVEC SEUIL LIMITE DE -1.5 AU MAXIMUM
comparison(df=Tab, velocity=Velocity, maptype='Mean Velocity', region='Main', plottype='MapDiff', type=0, vmin=-4, vmax=4, mapplot=0, thresh=-1.5)
# NUAGE DE POINTS DE DIFFERENCES SUR CARTE DE LOMBARDI
comparison(df=Tab, velocity=Velocity, maptype='Mean Velocity', region='Main', plottype='MapDiff', type=0, vmin=-4, vmax=4, mapplot='Lombardi', thresh=0)
# NUAGE DE POINTS DE DIFFERENCES EN VALEURS ABSOLUES SUR CARTE DE LOMBARDI
comparison(df=Tab, velocity=Velocity, maptype='Mean Velocity', region='Main', plottype='MapDiff', type='abs', vmin=0, vmax=4, mapplot='Lombardi', thresh=0)
# NUAGE DE POINTS DE DIFFERENCES AVEC SEUIL LIMITE DE 2 AU MINIMUM SUR CARTE DE LOMBARDI
comparison(df=Tab, velocity=Velocity, maptype='Mean Velocity', region='Main', plottype='MapDiff', type=0, vmin=-4, vmax=4, mapplot='Lombardi', thresh=2)
# NUAGE DE POINTS DE DIFFERENCES AVEC SEUIL LIMITE DE -2 AU MAXIMUM SUR CARTE DE LOMBARDI
comparison(df=Tab, velocity=Velocity, maptype='Mean Velocity', region='Main', plottype='MapDiff', type=0, vmin=-4, vmax=4, mapplot='Lombardi', thresh=-2)


## AFFICHAGE DES GRAPHIQUES DE DIFFERENCES
# NORMAL
comparison(df=Tab, velocity=Velocity, maptype='Mean Velocity', region='Main', plottype='Graph', type=0, vmin=0, vmax=0, mapplot=0, thresh=0)
# EN VALEUR ABSOLUE
comparison(df=Tab, velocity=Velocity, maptype='Mean Velocity', region='Main', plottype='Graph', type='abs', vmin=0, vmax=0, mapplot=0, thresh=0)


## AFFICHAGE DES CORRELATIONS
# NORMAL
comparison(df=Tab, velocity=Velocity, maptype='Mean Velocity', region='Main', plottype='Correlation', type=0, vmin=0, vmax=0, mapplot=0, thresh=0)
# EN VALEUR ABSOLUE
comparison(df=Tab, velocity=Velocity, maptype='Mean Velocity', region='Main', plottype='Correlation', type='abs', vmin=0, vmax=0, mapplot=0, thresh=0)


#=======================================================================================================================


### EXEMPLE 2 : CATALOGUE POUR ORION COMPLET
## CREATION DU CATALOGUE ET FILTRAGE
param_cornu = ['Pred', 'P(CI)', 'P(CII)']
param_kounkel = ['2MASS', 'Gaia', 'RVmean', 'e_RVmean', 'Av', 'e_Av']
param_johnsson = ['ID', 'HRV', 'errHRV']
param_gaia = ['source_id', 'parallax', 'parallax_error', 'pmra', 'pmra_error', 'pmdec', 'pmdec_error']
ProbThresh = 0
Apogee, Velocity, VelErr = 'Kounkel', 'RVmean', 'e_RVmean'
Remove = 'Both'
Loc = [0,  360,  -90,  90]

Tab = catalog(Param_cornu=param_cornu, Param_Gaia=param_gaia, Param_Johnsson=param_johnsson, Param_Kounkel=param_kounkel, ProbThresh=ProbThresh, Loc=Loc, RemoveType=Remove, Apo=Apogee)

Tab = position_filter(df=Tab, MinDist=300, MaxDist=500, error=50)
Tab = velocity_filter(df=Tab, velocity=Velocity, velerr=VelErr, vmin=1, vmax=16, error=0.5)
Tab = pm_filter(df=Tab,  pmramin=-10,  pmramax=10,  pmdecmin=-10,  pmdecmax=10,  type='relative',  error=0.5)


## AFFICHAGES 2D DES VITESSES RADIALES
# Remarque : Il faut bien modifier vmin et vmax pour fixer l'echelle de couleur qui dépend de l'échantillon de YSOs
# Remarque : Pour les deux fits, bien stipuler vinf et vsup pour séparer les deux échantillons
# Remarque : Il y a des problèmes que je n'ai pas pu résoudre à temps sur le nuage 3D de vecteurs autre que L1630S


# AFFICHAGE NUAGE DE POINTS SIMPLE
velocity_2D(df=Tab, velocity=Velocity, velerr=VelErr, pm=0, fit=0, order=0, vinf=0, vsup=0, vmin=7, vmax=13)
# AFFICHAGE NUAGE DE POINTS AVEC DEUX FITS D'ORDRE 1
velocity_2D(df=Tab, velocity=Velocity, velerr=VelErr, pm=0, fit=2, order=1, vinf=10, vsup=11, vmin=7, vmax=13)
# AFFICHAGE NUAGE DE POINTS AVEC DEUX FITS D'ORDRE 2
velocity_2D(df=Tab, velocity=Velocity, velerr=VelErr, pm=0, fit=2, order=2, vinf=10, vsup=11, vmin=7, vmax=13)
# AFFICHAGE NUAGE DE VECTEURS SIMPLE
velocity_2D(df=Tab, velocity=Velocity, velerr=VelErr, pm=1, fit=0, order=0, vinf=0, vsup=0, vmin=7, vmax=16)
# AFFICHAGE NUAGE DE VECTEURS AVEC DEUX FITS D'ORDRE 1
velocity_2D(df=Tab, velocity=Velocity, velerr=VelErr, pm=1, fit=2, order=1, vinf=10, vsup=11, vmin=7, vmax=13)
# AFFICHAGE NUAGE DE VECTEURS AVEC DEUX FITS D'ORDRE 2
velocity_2D(df=Tab, velocity=Velocity, velerr=VelErr, pm=1, fit=2, order=2, vinf=10, vsup=11, vmin=7, vmax=13)


## AFFICHAGES 3D DES VITESSES
# Remarque : Il faut bien modifier vmin et vmax pour fixer l'echelle de couleur qui dépend de l'échantillon de YSOs
# Remarque : Pour les deux fits, bien stipuler vinf et vsup pour séparer les deux échantillons
# Remarque : Pour le nuage de vecteurs, l'échelle de couleur concerne ici la vitesse totale (cf velocity_graphs)
# NUAGE DE POINTS SIMPLE
velocity_3D(df=Tab, velocity=Velocity, velerr=VelErr, pm=0, fit=0, order=0, vinf=0, vsup=0, vmin=7, vmax=13)
# NUAGE DE POINTS AVEC UN FIT D'ORDRE 1
velocity_3D(df=Tab, velocity=Velocity, velerr=VelErr, pm=0, fit=1, order=1, vinf=0, vsup=0, vmin=7, vmax=13)
# NUAGE DE POINTS AVEC UN FIT D'ORDRE 2
velocity_3D(df=Tab, velocity=Velocity, velerr=VelErr, pm=0, fit=1, order=2, vinf=0, vsup=0, vmin=7, vmax=13)
# NUAGE DE POINTS AVEC DEUX FITS D'ORDRE 1
velocity_3D(df=Tab, velocity=Velocity, velerr=VelErr, pm=0, fit=2, order=1, vinf=10, vsup=11, vmin=7, vmax=13)
# NUAGE DE POINTS AVEC DEUX FITS D'ORDRE 2
velocity_3D(df=Tab, velocity=Velocity, velerr=VelErr, pm=0, fit=2, order=2, vinf=10, vsup=11, vmin=7, vmax=13)


## EXEMPLES CLASSIQUES DE CARTOGRAPHIES
# NUAGE DE POINTS DE VITESSES RADIALES DANS ORION COMPLET
mapsplot(df=Tab, plottype='Scatter', maptype='Mean Velocity', region='Orion', p='Radial Velocity', velocity=Velocity, zoom=0, axis=0, vmin=1, vmax=16)
# NUAGE DE POINTS DE DISTANCES DANS ORION COMPLET
mapsplot(df=Tab, plottype='Scatter', maptype='Mean Velocity', region='Orion', p='Distance', velocity=Velocity, zoom=0, axis=0, vmin=320, vmax=450)
# NUAGE DE POINTS D'EXTINCTION DANS ORION COMPLET
mapsplot(df=Tab, plottype='Scatter', maptype='Mean Velocity', region='Orion', p='Extinction', velocity=Velocity, zoom=0, axis=0, vmin=10**20, vmax=10**22)
# NUAGE DE VECTEURS DE VITESSES RADIALES DANS ORION COMPLET
mapsplot(df=Tab, plottype='Quiver', maptype='Mean Velocity', region='Orion', p='Radial Velocity', velocity=Velocity, zoom=0, axis=0, vmin=1, vmax=16)
# NUAGE DE VECTEURS DE DISTANCES DANS ORION COMPLET
mapsplot(df=Tab, plottype='Quiver', maptype='Mean Velocity', region='Orion', p='Distance', velocity=Velocity, zoom=0, axis=0, vmin=320, vmax=450)
# NUAGE DE VECTEURS D'EXTINCTION DANS ORION COMPLET
mapsplot(df=Tab, plottype='Quiver', maptype='Mean Velocity', region='Orion', p='Extinction', velocity=Velocity, zoom=0, axis=0, vmin=10**20, vmax=10**22)


## AFFICHAGE DE CARTES DE DIFFERENCS DE VITESSES
# NUAGE DE POINTS DE DIFFERENCES
comparison(df=Tab, velocity=Velocity, maptype='Mean Velocity', region='Orion', plottype='MapDiff', type=0, vmin=-8, vmax=8, mapplot=0, thresh=0)
# NUAGE DE POINTS DE DIFFERENCES EN VALEURS ABSOLUES
comparison(df=Tab, velocity=Velocity, maptype='Mean Velocity', region='Orion', plottype='MapDiff', type='abs', vmin=0, vmax=8, mapplot=0, thresh=0)
# NUAGE DE POINTS DE DIFFERENCES AVEC SEUIL LIMITE DE 2 AU MINIMUM
comparison(df=Tab, velocity=Velocity, maptype='Mean Velocity', region='Orion', plottype='MapDiff', type=0, vmin=-8, vmax=8, mapplot=0, thresh=2)
# NUAGE DE POINTS DE DIFFERENCES AVEC SEUIL LIMITE DE -2 AU MAXIMUM
comparison(df=Tab, velocity=Velocity, maptype='Mean Velocity', region='Orion', plottype='MapDiff', type=0, vmin=-8, vmax=8, mapplot=0, thresh=-2)


## AFFICHAGE DES GRAPHIQUES DE DIFFERENCES
# NORMAL
comparison(df=Tab, velocity=Velocity, maptype='Mean Velocity', region='Orion', plottype='Graph', type=0, vmin=0, vmax=0, mapplot=0, thresh=0)
# EN VALEUR ABSOLUE
comparison(df=Tab, velocity=Velocity, maptype='Mean Velocity', region='Orion', plottype='Graph', type='abs', vmin=0, vmax=0, mapplot=0, thresh=0)


## AFFICHAGE DES CORRELATIONS
# NORMAL
comparison(df=Tab, velocity=Velocity, maptype='Mean Velocity', region='Orion', plottype='Correlation', type=0, vmin=0, vmax=0, mapplot=0, thresh=0)
# EN VALEUR ABSOLUE
comparison(df=Tab, velocity=Velocity, maptype='Mean Velocity', region='Orion', plottype='Correlation', type='abs', vmin=0, vmax=0, mapplot=0, thresh=0)
