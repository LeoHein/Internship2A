from Velocity import *

#=======================================================================================================================

# Retire tous les éléments d'une liste ayant une valeur donnée en entrée
## list : liste à modifier
## val : valeur à supprimer
### new_list : liste modifiée

def remove_values_from_list(list, val):
    new_list = [value for value in list if value != val]
    return new_list

#=======================================================================================================================

# Filtre les objets du Dataframe d'entrée, ne garde que les candidats ayant une probabilité supérieure au seuil donné
## df : Dataframe de travail
## ProbThresh : 0 pour utiliser les prédictions de Cornu ou valeur de probabilité pour le seuil minimal accepté
### df : Dataframe modifié

def Cornu_filter(df, ProbThresh):
    df = df[df['Pred'].isin([0, 1])]
    if ProbThresh != 0:
        df = df[(df['P(CI)'] > ProbThresh) | (df['P(CII)'] > ProbThresh)]
    return df

#=======================================================================================================================

# Filtre les objets du Dataframe d'entrée par leur position
## df  : Dataframe de travail
## Loc : liste en coordonnées équatoriale des 4 coordonnées limites [RA min, RA max, DEC min, DEC max] en degré
### df : Dataframe modifié

def Loc_filter(df, Loc):
    df = df[(df['RA'] > Loc[0]) & (df['RA'] < Loc[1])]
    df = df[(df['DEC'] > Loc[2]) & (df['DEC'] < Loc[3])]
    return df

#=======================================================================================================================

# Assigne à chaque objets du Dataframe de travail sa région correspondante parmi Orion-A, L1630S et L1630N
## df : Dataframe de travail
### df : Dataframe modifié

def assignation(df):
    df.loc[(df['RA'] > 86) & (df['RA'] < 88) & (df['DEC'] > -1) & (df['DEC'] < 1), 'Region'] = 'L1630N'
    df.loc[(df['RA'] > 85) & (df['RA'] < 86) & (df['DEC'] > -3) & (df['DEC'] < -0.8), 'Region'] = 'L1630S'
    df.loc[(df['RA'] > 83) & (df['RA'] < 86.5) & (df['DEC'] > -11) & (df['DEC'] < -4), 'Region'] = 'OrionA'
    return df

#=======================================================================================================================

# Permet de gérer le cas des étoiles doubles ou ayant changé de nomenclature entre DR2 et EDR3
## df : Dataframe de travail
## Param_Kounkel : paramètres du catalogue Kounkel que l'on souhaite garder dans le catalogue final
## Param_Gaia : paramètres du catalogue Gaia que l'on souhaite garder dans le catalogue final
## RemoveType : Type d'étoile à retirer, "NomencChange", "Double" ou "Both"
### tab : Dataframe modifié

def doublestars(df, Param_Kounkel, Param_Gaia, RemoveType):
    dfc = df.copy()
    dfc = dfc.replace(r'^\s*$', np.NaN, regex=True)
    dfc.reset_index(drop=True, inplace=True)
    dfc.drop(dfc[(dfc['Gaia'].notnull()) & (dfc['source_id'].notnull())].index, inplace=True)
    dg = blank_cleaning(df, ['source_id', 'Gaia'])
    dg['erreur'] = dg['Gaia'] - dg['source_id']
    dg.loc[(dg['Gaia'] - dg['source_id']).abs() != 0, 'Type'] = 'Double'
    dg.loc[((dg['Gaia'] - dg['source_id']).abs() < 4310000000) & ((dg['Gaia'] - dg['source_id']).abs() > 4290000000),
           'Type'] = 'NomencChange'

    if RemoveType == 'Both':
        dg = dg[dg['Type'] != 'Double']
        dg = dg[dg['Type'] != 'NomChange']

    if RemoveType == 'NomencChange':
        dg = dg[dg['Type'] != 'NomChange']

    if RemoveType == 'double':
        dg = dg[dg['Type'] != 'Double']

    if RemoveType == 'None':
        mask = dg['Type'] == 'Double'
        dg = dg.reindex(dg[mask].index.repeat(2)).append(dg[~mask])
        L1 = dg.duplicated(keep='first')
        L2 = dg.duplicated(keep='last')
        Param_Kounkel.remove('RA')
        Param_Kounkel.remove('DEC')
        Param_Gaia.remove('RA')
        Param_Gaia.remove('DEC')
        dg.loc[L1, Param_Kounkel] = np.nan
        dg.loc[L2, Param_Gaia] = np.nan

    dfc.reset_index(drop=True, inplace=True)
    dg.reset_index(drop=True, inplace=True)
    tab = pd.concat([dfc, dg], axis=0, ignore_index=True)
    return tab

#=======================================================================================================================

# Créé un Dataframe contenant les données du catalogue de Cornu avec filtrage sur plusieurs variables
## ProbThresh : seuil de probabilité minimal accepté
## Loc : liste en coordonnées équatoriale des 4 coordonnées limites [RA min, RA max, DEC min, DEC max] en degré
## param : paramètres du catalogue de Cornu que l'on souhaite garder dans le catalogue final
### cornu : nouveau catalogue de Cornu

def import_cornu(ProbThresh, Loc, param):
    cornu = pd.read_csv('Cornu.txt', engine='python', sep='|', skiprows=10, skipfooter=1, header=None, names=(
                        'RAdegDEdeg', 'catalog', 'OriClass', '3.6mag', 'e_3.6mag', '4.5mag', 'e_4.5mag', '5.8mag',
                        'e_5.8mag', '8.0mag', 'e_8.0mag', '24mag', 'e_24mag', 'target', 'Pred', 'P(CI)', 'P(CII)',
                        'P(Other)'))
    cornu.drop_duplicates(subset="RAdegDEdeg", keep='first', inplace=True)
    cornu[['RA', 'DEC']] = cornu.RAdegDEdeg.str.split(expand=True)
    cornu[['RA', 'DEC']] = cornu[['RA', 'DEC']].apply(pd.to_numeric)
    del cornu['RAdegDEdeg']
    cornu = Cornu_filter(cornu, ProbThresh)
    cornu = Loc_filter(cornu, Loc)
    cornu = cornu[param]
    cornu.to_csv('cornu_filtered.csv', header=True, decimal='.', index=False, sep='\t', float_format='%.7f')
    cornu.reset_index(drop=True, inplace=True)
    return cornu

#=======================================================================================================================

# Créé un Dataframe contenant les données du catalogue de Kounkel avec filtrage sur les paramètres
## param : paramètres du catalogue de Kounkel que l'on souhaite garder dans le catalogue final
### Kounkel : nouveau catalogue de Kounkel

def import_kounkel(param):
    Kounkel = pd.read_csv('Kounkel.csv', engine='python', sep=';', skiprows=0, skipfooter=0, header=0)
    Kounkel = Kounkel.rename(columns={'RAJ2000': 'RA', 'DEJ2000': 'DEC'})
    Kounkel.drop_duplicates(subset="2MASS", keep='first', inplace=True)
    Kounkel = Kounkel[param]
    return Kounkel

#=======================================================================================================================

# Créé un Dataframe contenant les données du crossmatch Cornu x Johnsson avec filtrage selon plusieurs variables
## ProbThresh : seuil de probabilité minimal accepté
## Loc : liste en coordonnées équatoriale des 4 coordonnées limites [RA min, RA max, DEC min, DEC max] en degré
## param : paramètres du catalogue de Cornu et de Johnsson que l'on souhaite garder dans le catalogue final
### cornu : catalogue crossmatch Cornu X Johnsson

def import_johnsson(ProbThresh, Loc, param):
    cornuXApodr16 = pd.read_csv('cornuXApodr16.csv', engine='python', sep=',', skiprows=0, skipfooter=0, header=0)
    cornuXApodr16 = cornuXApodr16.rename(columns={'RAJ2000': 'RA', 'DEJ2000': 'DEC'})
    cornuXApodr16.drop_duplicates(subset="ID", keep='first', inplace=True)
    cornuXApodr16 = Cornu_filter(cornuXApodr16, ProbThresh)
    cornuXApodr16 = Loc_filter(cornuXApodr16, Loc)
    cornuXApodr16 = cornuXApodr16[param]
    cornuXApodr16.reset_index(drop=True, inplace=True)
    return cornuXApodr16

#=======================================================================================================================

# Créé un Dataframe contenant les données du crossmatch entre deux Dataframes d'entrée
## df1 : Dataframe 1 à crossmatcher
## df2 : Dataframe 2 à crossmatcher
## reste : 0 pour crossmatch exclusif et 1 pour crossmatch inclusif
### tab : Dataframe contenant les données crossmatchées des deux paramètres d'entrée

def crossmatch(df1, df2, reste):
    dftemp = df1.copy()
    par1 = df1.columns
    par2 = df2.columns
    from astroML.crossmatch import crossmatch_angular
    coords1_equ = df1[['RA', 'DEC']].to_numpy()
    coords2_equ = df2[['RA', 'DEC']].to_numpy()

    L1, L2, L3, L4 = [], [], [], []
    for c in range(len(coords1_equ)):
        L1.append(coords1_equ[c][0])
        L2.append(coords1_equ[c][1])
    for d in range(len(coords2_equ)):
        L3.append(coords2_equ[d][0])
        L4.append(coords2_equ[d][1])
    coords1_equ = [L1, L2]
    coords2_equ = [L3, L4]

    imX = np.zeros((len(coords1_equ[0]), 2))
    stX = np.zeros((len(coords2_equ[0]), 2))
    imX[:, 0] = coords1_equ[0]
    imX[:, 1] = coords1_equ[1]
    stX[:, 0] = coords2_equ[0]
    stX[:, 1] = coords2_equ[1]
    max_radius = 1 / 3600  # 1 arcsec
    dist, ind = crossmatch_angular(imX, stX, max_radius)
    pas_trouves = np.where(ind == len(coords2_equ[0]))
    df1.reset_index(drop=True, inplace=True)
    df1.drop(pas_trouves[0].tolist(), 0, inplace=True)
    ind = ind.tolist()
    ind = remove_values_from_list(ind, len(coords2_equ[0]))
    tab = df2.iloc[ind]
    tab.reset_index(drop=True, inplace=True)
    df1.reset_index(drop=True, inplace=True)
    tab = pd.concat([tab, df1], axis=1, ignore_index=True)
    tab.columns = list(par2) + list(par1)
    tab = tab.loc[:, ~tab.columns.duplicated()]
    if reste == 1:
        dftemp = dftemp.iloc[pas_trouves]
        tab = pd.concat([tab, dftemp], axis=0)
        rows = df2.index[ind]
        df2.drop(rows, inplace=True)
        tab = pd.concat([tab, df2], axis=0)
    return tab

#=======================================================================================================================

# Fonction principale pour la création du catalogue
## Param_cornu : paramètres du catalogue de Cornu que l'on souhaite garder dans le catalogue final
## Param_Gaia : paramètres du catalogue de Gaia que l'on souhaite garder dans le catalogue final
## Param_Johnsson : paramètres du catalogue de Johnsson que l'on souhaite garder dans le catalogue final
## Param_Kounkel : paramètres du catalogue de Kounkel que l'on souhaite garder dans le catalogue final
## ProbThresh : seuil de probabilité minimal accepté
## Loc : liste en coordonnées équatoriale des 4 coordonnées limites [RA min, RA max, DEC min, DEC max] en degré
## RemoveType : Type d'étoile à retirer, "NomencChange", "Double" ou "Both"
## Apo : choix du catalogue APOGEE : "Kounkel", "Johnsson" ou "Both"
### CornuXapoXgaia : catalogue final

def catalog(Param_Cornu, Param_Gaia, Param_Johnsson, Param_Kounkel, ProbThresh, Loc, RemoveType, Apo):
    for L in [Param_Cornu, Param_Gaia, Param_Johnsson, Param_Kounkel]:
        L.append('RA')
        L.append('DEC')
    gaiak = True
    gaiag = True
    if 'Gaia' not in Param_Kounkel:
        gaiak = False
        Param_Kounkel.append('Gaia')
    if 'source_id' not in Param_Gaia:
        gaiag = False
        Param_Gaia.append('source_id')

    cornuXedr3 = pd.read_csv('cornuXedr3.csv', engine='python', sep=',', skiprows=0, skipfooter=0, header=0)
    cornuXedr3 = Cornu_filter(cornuXedr3, ProbThresh)
    cornuXedr3 = cornuXedr3.rename(columns={'RAJ2000': 'RA', 'DEJ2000': 'DEC'})
    pd.set_option('display.max_columns', None)
    cornuXedr3 = Loc_filter(cornuXedr3, Loc)
    cornuXedr3.drop_duplicates(subset="source_id", keep='first', inplace=True)
    cornuXedr3 = cornuXedr3.astype({"source_id": str})
    cornuXedr3 = cornuXedr3[list(set(Param_Gaia + Param_Cornu))]
    if 'parallax' in Param_Gaia:
        cornuXedr3['ProbDist'] = 1000 / cornuXedr3['parallax']
        Param_Gaia.append('ProbDist')
        if 'parallax_error' in Param_Gaia:
            cornuXedr3['MinDist'] = 1000 / (cornuXedr3['parallax'] + cornuXedr3['parallax_error'])
            cornuXedr3['MaxDist'] = 1000 / (cornuXedr3['parallax'] - cornuXedr3['parallax_error'])
            Param_Gaia.append('MaxDist')
            Param_Gaia.append('MinDist')
    cornuXedr3.to_csv('res_cornuXedr3.txt', header=True, decimal='.', index=False, sep='\t', float_format='%.7f')

    cornu = import_cornu(ProbThresh, Loc, Param_Cornu)
    if Apo == 'Kounkel':
        apogee = import_kounkel(Param_Kounkel)
        cornuXapo = crossmatch(cornu, apogee, 0)
    if Apo == 'Johnsson':
        cornuXapo = import_johnsson(ProbThresh, Loc, Param_Johnsson)
    if Apo == 'Both':
        apogee = import_kounkel(Param_Kounkel)
        cornuXkounkel = crossmatch(cornu, apogee, 0)
        cornuXjohnsson = import_johnsson(ProbThresh, Loc, Param_Johnsson)
        cornuXjohnsson.reset_index(drop=True, inplace=True)
        cornuXkounkel.reset_index(drop=True, inplace=True)
        cornuXapo = crossmatch(cornuXjohnsson, cornuXkounkel, 0)

    cornuXapo.to_csv('res_cornuXapo.txt', header=True, decimal='.', index=False, sep='\t', float_format='%.7f')
    cornuXapoXgaia = crossmatch(cornuXapo, cornuXedr3, 1)
    cornuXapoXgaia = assignation(cornuXapoXgaia)
    if Apo == 'Kounkel':
        cornuXapoXgaia = doublestars(cornuXapoXgaia, Param_Kounkel, Param_Gaia, RemoveType)
    if not gaiak:
        cornuXapoXgaia.drop(columns='Gaia')
    if not gaiag:
        cornuXapoXgaia.drop(columns='source_id')
    return cornuXapoXgaia



