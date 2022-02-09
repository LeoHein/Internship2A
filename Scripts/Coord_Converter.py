from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, LSR
import pandas as pd

#=======================================================================================================================

# Convertit des coordonnées équatoriales en coordonnées galactiques
## RA : ascension droite en degré
## DEC : déclinaison en degré
### glon : longitude galactique en degré
### glat : latitude galactique en degré

def EquToGal(RA, DEC):
    hcg7_center = SkyCoord(RA * u.degree, DEC * u.degree, frame='icrs')
    glon = hcg7_center.galactic.l.deg
    glat = hcg7_center.galactic.b.deg
    return glon, glat

#=======================================================================================================================

# Convertit des coordonnées galactiques en coordonnées équatoriales
## glon : longitude galactique en degré
## glat : latitude galactique en degré
### RA : ascension droite en degré
### DEC : déclinaison en degré

def GalToEqu(glon, glat):
    c = SkyCoord(l=glon * u.degree, b=glat * u.degree, frame='galactic')
    c_fk5 = c.transform_to('fk5')
    RA = c_fk5.ra.deg
    DEC = c_fk5.dec.deg
    return RA, DEC

#=======================================================================================================================

# Convertit une vitesse radiale d'un corps dans le référentiel héliocentrique au référentiel LSR
## RA : ascension droite de l'objet en degré
## DEC : déclinaison de l'objet en degré
## v : vitesse radiale de l'objet en km/s
### new_rv : vitesse radiale de l'objet dans le référentiel LSR en km/s

def to_LSR(RA, DEC, v):
    coords = [RA, DEC]
    my_observation = ICRS(ra=coords[0]*u.deg, dec=coords[1]*u.deg,
                          pm_ra_cosdec=0*u.mas/u.yr, pm_dec=0*u.mas/u.yr,
                          radial_velocity=v*u.km/u.s, distance=1*u.pc)
    new_rv = my_observation.transform_to(LSR()).radial_velocity
    return new_rv

#=======================================================================================================================

# Ajoute les coordonnées galactiques au dataframe d'entrée à partir des coordonnées équatoriales
## df : Dataframe dans lequel ajouter les colonnes ['glon'] et ['glat'] des coordonnées galactiques
## RA : Nom de la colonne de df contenant les ascensions droites
## DEC : Nom de la colonne de df contenant les déclinaisons
### df : Dataframe modifié

def Gal_implement(df, RA, DEC):
    df['gal'] = df.apply(lambda x: EquToGal(x[RA], x[DEC]), axis=1)
    df[['glon', 'glat']] = pd.DataFrame(df['gal'].tolist(), index=df.index)
    del df['gal']
    return df

