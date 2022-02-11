import matplotlib.pyplot as plt
import matplotlib.cm as cmx
from mpl_toolkits.mplot3d import *
import matplotlib
from Coord_Converter import *
from Fitting import *

#=======================================================================================================================

# Affiche le nuage de point des vitesses radiales en 2D avec plusieurs ajouts possibles
## df : Dataframe de travail
## velocity : nom du paramètre vitesse radiale
## velerr : nom du paramètre erreur de la vitesse radiale
## pm : 0 sans ou 1 avec les mouvements propres (nuage de vecteurs)
## fit : 0 sans, 2 avec deux surfaces (une pour chaque échantillon)
## order : ordre de la courbe de fitting : 1 pour linéaire et 2 pour quadratique
## vinf, vsup : respectivement les bornes supérieures et inférieures des deux échantillons pour le fit en km/s
## vmin, vmax : respectivement les vitesses minimales et maximales pour l'échelle de couleur en km/s

def velocity_2D(df, velocity, velerr, pm, fit, order, vinf, vsup, vmin, vmax):
    df = vLSR(df, velocity)
    velLSR = velocity+'LSR'
    df = blank_cleaning(df, [velLSR, velerr])
    if pm == 0:
        x = df['RA']
        y = df['DEC']
        v = df[velLSR]
        cm = plt.get_cmap('coolwarm')
        cNorm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
        plt.xlabel('RA (deg)')
        plt.ylabel('DEC (deg)')
        plt.scatter(x, y, c=scalarMap.to_rgba(v))
        ax = plt.gca()
        ax.set_xlim(ax.get_xlim()[::-1])
        scalarMap.set_array(v)
        cbar = plt.colorbar(scalarMap)
        cbar.set_label('YSOs Radial Velocity (km/s)', rotation=270, labelpad=15)
    if pm == 1:
        df = blank_cleaning(df, ['pmra', 'pmdec'])
        x = df['RA']
        y = df['DEC']
        colors = df[velLSR]
        u = -df['pmra']
        v = df['pmdec']
        cm = cmx.coolwarm
        cNorm = matplotlib.colors.Normalize(vmin, vmax)
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
        plt.xlabel('RA (deg)')
        plt.ylabel('DEC (deg)')
        plt.quiver(x, y, u, v, color=cm(cNorm(colors)), edgecolor='black', linewidth=0.2)
        ax = plt.gca()
        ax.set_xlim(ax.get_xlim()[::-1])
        scalarMap.set_array(v)
        cbar = plt.colorbar(scalarMap)
        cbar.set_label('YSOs Radial Velocity (km/s)', rotation=270, labelpad=15)

    if fit == 2:
        df1 = df[df[velLSR] < vinf]
        df2 = df[df[velLSR] > vsup]
        x1 = df1['RA']
        y1 = df1['DEC']
        x2 = df2['RA']
        y2 = df2['DEC']
        abs = np.linspace(min(x), max(x), 100)
        if order == 1:
            courbe = line_fitting(x1, y1, order)
            plt.plot(abs, courbe[0] + courbe[1] * abs, c=cm(cNorm(df1[velLSR].mean())))
            courbe = line_fitting(x2, y2, order)
            plt.plot(abs, courbe[0] + courbe[1] * abs, c=cm(cNorm(df2[velLSR].mean())))
        if order == 2:
            courbe = line_fitting(x1, y1, order)
            plt.plot(abs, courbe[0] + courbe[1] * abs + courbe[2] * abs ** 2, c=cm(cNorm(df1[velLSR].mean())))
            courbe = line_fitting(x2, y2, order)
            plt.plot(abs, courbe[0] + courbe[1] * abs + courbe[2] * abs ** 2, c=cm(cNorm(df2[velLSR].mean())))
    #plt.axis([min(x), max(x), min(y), max(y)])
    plt.show()

#=======================================================================================================================

# Affiche le nuage de point des vitesses radiales en 3D avec plusieurs ajouts possibles
## df : Dataframe de travail
## velocity : nom du paramètre vitesse radiale
## velerr : nom du paramètre erreur de la vitesse radiale
## pm : 0 sans ou 1 avec les mouvements propres (nuage de vecteurs)
## fit : 0 sans, 1 avec une surface moyenne totale, 2 avec deux surfaces (une pour chaque échantillon)
## order : ordre de la surface de fitting : 1 pour linéaire et 2 pour quadratique
## vinf, vsup : respectivement les bornes supérieures et inférieures des deux échantillons pour le fit en km/s
## vmin, vmax : respectivement les vitesses minimales et maximales pour l'échelle de couleur en km/s

def velocity_3D(df, velocity, velerr, pm, fit, order, vinf, vsup, vmin, vmax):
    df = vLSR(df, velocity)
    velLSR = velocity+'LSR'
    df = blank_cleaning(df, [velLSR, velerr, 'ProbDist'])
    if pm == 0:
        x = df['RA']
        y = df['DEC']
        z = df['ProbDist']
        v = df[velLSR]
        cm = plt.get_cmap('coolwarm')
        cNorm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.scatter(x, y, z, c=scalarMap.to_rgba(v), vmin=vmin, vmax=vmax)
        scalarMap.set_array(v)
        cbar = fig.colorbar(scalarMap, shrink=0.6)
        cbar.set_label('YSOs Radial Velocity (km/s)', rotation=270, labelpad=15)
        ax.set_xlabel('RA')
        ax.set_ylabel('DEC')
        ax.set_zlabel('ProbDist')
        ax = plt.gca()
        ax.set_xlim(ax.get_xlim()[::-1])
    if pm == 1:
        df = blank_cleaning(df, ['pmra', 'pmdec'])
        df['VRA'] = 3.08 * 10 ** 13 / (24 * 3600 * 365.25) * df['ProbDist'] * np.pi / (3600 * 1000 * 180) * df['pmra']
        df['VDEC'] = 3.08 * 10 ** 13 / (24 * 3600 * 365.25) * df['ProbDist'] * np.pi / (3600 * 1000 * 180) * df['pmdec']
        df = vLSR(df, 'VRA')
        df = vLSR(df, 'VDEC')
        df[velLSR] = df[velLSR] - df[velLSR].mean()
        df['VRALSR'] = df['VRALSR'] - df['VRALSR'].mean()
        df['VDECLSR'] = df['VDECLSR'] - df['VDECLSR'].mean()
        df['VTOT'] = np.sqrt(df['VRALSR'] ** 2 + df['VDECLSR'] ** 2 + df[velLSR] ** 2)

        x = df['RA']
        y = df['DEC']
        z = df['ProbDist']
        u = df['VRALSR']*0.01
        v = df['VDECLSR']*0.02
        w = df[velLSR]*10
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        cm = plt.get_cmap('coolwarm')
        cNorm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
        ax.scatter(x, y, z, c=scalarMap.to_rgba(df['VTOT']), vmin=vmin, vmax=vmax, s=10)
        ax.set_xlabel('RA')
        ax.set_ylabel('DEC')
        ax.set_zlabel('ProbDist')
        cm = cmx.coolwarm
        colors = df['VTOT']
        cNorm = matplotlib.colors.Normalize(vmin, vmax)
        ax.quiver(x, y, z, u, v, w, color=cm(cNorm(colors)), arrow_length_ratio=0.00000005, pivot='tail')
        cbar = fig.colorbar(scalarMap, shrink=0.6)
        cbar.set_label('YSOs Radial Velocity (km/s)', rotation=270, labelpad=15)
        ax = plt.gca()
        ax.set_xlim(ax.get_xlim()[::-1])
    if fit == 1:
        X, Y, Z = plane_fitting(x, y, z, order)
        ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.2)
    if fit == 2:
        df1 = df[df[velLSR] < vinf]
        df2 = df[df[velLSR] > vsup]
        X, Y, Z = plane_fitting(df1['RA'], df1['DEC'], df1['ProbDist'], order)
        ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.2, color='blue')
        U, V, W = plane_fitting(df2['RA'], df2['DEC'], df2['ProbDist'], order)
        ax.plot_surface(U, V, W, rstride=1, cstride=1, alpha=0.2, color='red')
        print(df1.shape[0], "YSO avec vitesse radiale inférieure à", vinf, "km/s")
        print(df2.shape[0], "YSO avec vitesse radiale supérieure à", vsup, "km/s")
    plt.show()

#=======================================================================================================================

# Ajoute au Dataframe de travail une colonne en transformant les vitesses radaiales dans le système LSR
## df : Dataframe de travail
## velocity : nom du paramètre vitesse radiale
### df : Dataframe modifié

def vLSR(df, velocity):
    coords = df[['RA', 'DEC', velocity]].to_numpy()
    velLSR = to_LSR(coords[:, 0], coords[:, 1], coords[:, 2])
    name = velocity + 'LSR'
    df[name] = velLSR
    return df

#=======================================================================================================================

# Supprime du Dataframe de travail les objets pour lesquels il manque les données d'au moins un paramètre de la liste L
## df : Dataframe de travail
## L : liste des paramètres à nettoyer
### df : Dataframe modifié

def blank_cleaning(df, L):
    for param in L:
        df[param] = pd.to_numeric(df[param], errors='coerce')
        df.dropna(subset=[param], inplace=True)
    return df

#=======================================================================================================================

# Filtre les objets du Dataframe d'entrée selon leurs vitesses radiales et leurs erreurs associées
## df : Dataframe de travail
## velocity : nom du paramètre vitesse radiale
## velerr : nom du paramètre erreur de la vitesse radiale
## vmin, vmax : respectivement bornes inférieure et supérieure de l'intervalle de vitesses radiales accepté en km/s
## error : 0 pour aucun filtre sur l'erreur ou un nombre positif comme seuil d'erreur maximum accepté en km/s
### df : Dataframe modifié

def velocity_filter(df, velocity, velerr, vmin, vmax, error):
    df = vLSR(df, velocity)
    velLSR = velocity + 'LSR'
    df = blank_cleaning(df, [velLSR, velerr])
    if error == 0:
        df = df[(df[velLSR] > vmin) & (df[velLSR] < vmax)]
    else:
        df = df[(df[velLSR] > vmin) & (df[velLSR] < vmax)]
        df = df[df[velerr] < error]
    return df

#=======================================================================================================================

# Filtre les objets du Dataframe d'entrée selon leurs mouvements propres et leurs erreurs associées
## df : Dataframe de travail
## pmramin, pmramax : respectivement bornes inférieure et supérieure de l'intervalle de mp selon ra en mas/year
## pmdecmin, pmdecmax : respectivement bornes inférieure et supérieure de l'intervalle de mp selon dec en mas/year
## type : type d'erreur, "uniform" pour une érreur linéaire, "relative" pour relative
## error : 0 pour aucun filtre sur l'erreur ou un nombre positif comme seuil d'erreur maximum accepté en mas/year
### df : Dataframe modifié

def pm_filter(df, pmramin, pmramax, pmdecmin, pmdecmax, type, error):
    df = blank_cleaning(df, ['pmra', 'pmdec', 'pmra_error', 'pmdec_error'])
    if type == 'uniform':
        if error == 0:
            df = df[(df['pmra'] > pmramin) & (df['pmra'] < pmramax)]
            df = df[(df['pmdec'] > pmdecmin) & (df['pmra'] < pmdecmax)]
        else:
            df = df[(df['pmra'] > pmramin) & (df['pmra'] < pmramax)]
            df = df[df['pmra_error'] < error]
            df = df[(df['pmdec'] > pmdecmin) & (df['pmdec'] < pmdecmax)]
            df = df[df['pmdec_error'] < error]

    if type == 'relative':
        if error == 0:
            df = df[(df['pmra'] > pmramin) & (df['pmra'] < pmramax)]
            df = df[(df['pmdec'] > pmdecmin) & (df['pmra'] < pmdecmax)]
        else:
            df = df[(df['pmra'] > pmramin) & (df['pmra'] < pmramax)]
            df = df[abs(df['pmra_error']) / abs(df['pmra']) < error]
            df = df[(df['pmdec'] > pmdecmin) & (df['pmdec'] < pmdecmax)]
            df = df[abs(df['pmdec_error']) / abs(df['pmdec']) < error]
    return df

#=======================================================================================================================

# Affiche par ordre croissant les vitesses radiales et erreurs associées des objets dans le Dataframe d'entrée
## df : Dataframe de travail
## velocity : nom du paramètre vitesse radiale
## velerr : nom du paramètre erreur de la vitesse radiale
## type : "Relative" pour enlever le mouvement radial moyen du nuage

def velocity_graph(df, velocity, velerr, type):
    df = vLSR(df, velocity)
    velLSR = velocity + 'LSR'
    df = blank_cleaning(df, [velLSR, velerr])
    df = df.sort_values(by=velLSR, ascending=True)
    y = df[velLSR]
    if type == 'Relative':
        y = df[velLSR] - df[velLSR].mean()
    x = [i for i in range(len(df))]
    err = df[velerr]
    plt.errorbar(x, y, yerr=err, fmt='.', color='black', ecolor='lightgray', elinewidth=3, capsize=0, markersize=2)
    plt.plot(np.linspace(0, df.shape[0], 2), [df[velLSR].mean()] * 2, label='Mean = ' + str(round(df[velLSR].mean(), 3)) + ' km/s')
    plt.plot(np.linspace(0, df.shape[0], 2), [df[velLSR].quantile(q=0.25)] * 2, label='1st quartile = ' + str(round(df[velLSR].quantile(0.25), 3)) + ' km/s')
    plt.plot(np.linspace(0, df.shape[0], 2), [df[velLSR].quantile(q=0.75)] * 2, label='3rd quartile = ' + str(round(df[velLSR].quantile(0.75), 3)) + ' km/s')
    plt.xlabel("YSO's ID")
    plt.ylabel("Radial velocity + error (km/s)")
    plt.title("YSO's radial velocity in ascending order")
    plt.legend()
    plt.show()

#=======================================================================================================================

# Affiche par ordre croissant les mouvements propres et erreurs associées des objets dans le Dataframe d'entrée
## df : Dataframe de travail
## abs : 1 pour afficher les données en valeurs absolues

def pm_graph(df, abs):
    df = blank_cleaning(df, ['pmra', 'pmdec', 'pmra_error', 'pmdec_error'])
    fig, axs = plt.subplots(1, 2, sharey=True)
    x = [i for i in range(len(df))]

    df = df.sort_values(by='pmra', ascending=True)
    if abs == 1:
        df['pmra'] = df['pmra'].abs()
        df = df.sort_values(by='pmra', ascending=True)
    yra = df['pmra']
    err = df['pmra_error']
    axs[0].errorbar(x, yra, yerr=err, fmt='.', color='black', ecolor='lightgray', elinewidth=3, capsize=0, markersize=2)
    axs[0].plot(np.linspace(0, df.shape[0], 2), [df['pmra'].mean()] * 2, label='Mean = ' + str(round(df['pmra'].mean(), 3)) + ' km/s')
    axs[0].plot(np.linspace(0, df.shape[0], 2), [df['pmra'].quantile(q=0.25)] * 2, label='1st quartile = ' + str(round(df['pmra'].quantile(0.25), 3)) + ' km/s')
    axs[0].plot(np.linspace(0, df.shape[0], 2), [df['pmra'].quantile(q=0.75)] * 2, label='3rd quartile = ' + str(round(df['pmra'].quantile(0.75), 3)) + ' km/s')
    axs[0].legend()

    df = df.sort_values(by='pmdec', ascending=True)
    if abs == 1:
        df['pmdec'] = df['pmdec'].abs()
        df = df.sort_values(by='pmdec', ascending=True)
    ydec = df['pmdec']
    err = df['pmdec_error']
    axs[1].errorbar(x, ydec, yerr=err, fmt='.', color='black', ecolor='lightgray', elinewidth=3, capsize=0, markersize=2)
    axs[1].plot(np.linspace(0, df.shape[0], 2), [df['pmdec'].mean()] * 2, label='Mean = ' + str(round(df['pmdec'].mean(), 3)) + ' km/s')
    axs[1].plot(np.linspace(0, df.shape[0], 2), [df['pmdec'].quantile(q=0.25)] * 2, label='1st quartile = ' + str(round(df['pmdec'].quantile(0.25), 3)) + ' km/s')
    axs[1].plot(np.linspace(0, df.shape[0], 2), [df['pmdec'].quantile(q=0.75)] * 2, label='3rd quartile = ' + str(round(df['pmdec'].quantile(0.75), 3)) + ' km/s')
    fig.suptitle("YSO's proper motions in ascending order (RA & DEC)")
    axs[1].legend()
    plt.show()

#=======================================================================================================================

# Affiche par ordre croissant toutes les vitesses et erreurs associées des objets dans le Dataframe d'entrée
## df : Dataframe de travail
## velocity : nom du paramètre vitesse radiale
## velerr : nom du paramètre erreur de la vitesse radiale
## type : 'Relative' pour retirer les mouvements moyens du nuage

def velocity_graphs(df, velocity, velerr, type):
    df = vLSR(df, velocity)
    velLSR = velocity +'LSR'
    df = blank_cleaning(df, [velLSR, velerr, 'pmra', 'pmdec', 'pmra_error', 'pmdec_error', 'ProbDist'])
    df['VRA'] = 3.08 * 10 ** 13 / (24 * 3600 * 365.25) * df['ProbDist'] * np.pi / (3600 * 1000 * 180) * df['pmra']
    df['eVRA'] = 3.08 * 10 ** 13 / (24 * 3600 * 365.25) * df['ProbDist'] * np.pi / (3600 * 1000 * 180) * df['pmra_error']
    df['VDEC'] = 3.08 * 10 ** 13 / (24 * 3600 * 365.25) * df['ProbDist'] * np.pi / (3600 * 1000 * 180) * df['pmdec']
    df['eVDEC'] = 3.08 * 10 ** 13 / (24 * 3600 * 365.25) * df['ProbDist'] * np.pi / (3600 * 1000 * 180) * df['pmdec_error']
    df = vLSR(df, 'VRA')
    df = vLSR(df, 'VDEC')
    df['VTOT'] = np.sqrt(df['VRALSR'] ** 2 + df['VDECLSR'] ** 2 + df[velLSR] ** 2)
    df['eVTOT'] = np.sqrt(df['eVRA'] ** 2 + df['eVDEC'] ** 2 + df[velerr] ** 2)
    fig, axs = plt.subplots(2, 2)
    x = [i for i in range(len(df))]

    if type == 'Relative':
        df[velLSR] = df[velLSR] - df[velLSR].mean()
        df['VRALSR'] = df['VRALSR'] - df['VRALSR'].mean()
        df['VDECLSR'] = df['VDECLSR'] - df['VDECLSR'].mean()
        df['VTOT'] = np.sqrt(df['VRALSR'] ** 2 + df['VDECLSR'] ** 2 + df[velLSR] ** 2)

    yvra = df['VRALSR'].sort_values(ascending=True)
    err = df['eVRA']
    axs[0, 0].errorbar(x, yvra, yerr=err, fmt='.', color='black', ecolor='lightgray', elinewidth=3, capsize=0, markersize=2)
    axs[0, 0].set_title("RA velocity")
    axs[0, 0].set_xlabel("YSOs ID")
    axs[0, 0].set_ylabel("Velocity (km/s)")
    axs[0, 0].plot(np.linspace(0, df.shape[0], 2), [df['VRALSR'].mean()] * 2)
    axs[0, 0].plot(np.linspace(0, df.shape[0], 2), [df['VRALSR'].mean()] * 2,
             label='Mean = ' + str(round(df[velLSR].mean(), 3)) + ' km/s')
    axs[0, 0].plot(np.linspace(0, df.shape[0], 2), [df['VRALSR'].quantile(q=0.25)] * 2,
             label='1st quartile = ' + str(round(df['VRALSR'].quantile(0.25), 3)) + ' km/s')
    axs[0, 0].plot(np.linspace(0, df.shape[0], 2), [df['VRALSR'].quantile(q=0.75)] * 2,
             label='3rd quartile = ' + str(round(df['VRALSR'].quantile(0.75), 3)) + ' km/s')
    axs[0, 0].legend()

    ydec = df['VDECLSR'].sort_values(ascending=True)
    err = df['eVDEC']
    axs[0, 1].errorbar(x, ydec, yerr=err, fmt='.', color='black', ecolor='lightgray', elinewidth=3, capsize=0, markersize=2)
    axs[0, 1].set_title("DEC velocity")
    axs[0, 1].set_xlabel("YSOs ID")
    axs[0, 1].set_ylabel("Velocity (km/s)")
    axs[0, 1].plot(np.linspace(0, df.shape[0], 2), [df['VDECLSR'].mean()] * 2)
    axs[0, 1].plot(np.linspace(0, df.shape[0], 2), [df['VDECLSR'].mean()] * 2,
             label='Mean = ' + str(round(df['VDECLSR'].mean(), 3)) + ' km/s')
    axs[0, 1].plot(np.linspace(0, df.shape[0], 2), [df['VDECLSR'].quantile(q=0.25)] * 2,
             label='1st quartile = ' + str(round(df['VDECLSR'].quantile(0.25), 3)) + ' km/s')
    axs[0, 1].plot(np.linspace(0, df.shape[0], 2), [df['VDECLSR'].quantile(q=0.75)] * 2,
             label='3rd quartile = ' + str(round(df['VDECLSR'].quantile(0.75), 3)) + ' km/s')
    axs[0, 1].legend()

    yrv = df[velLSR].sort_values(ascending=True)
    err = df[velerr]
    axs[1, 0].errorbar(x, yrv, yerr=err, fmt='.', color='black', ecolor='lightgray', elinewidth=3, capsize=0, markersize=2)
    axs[1, 0].set_title("Radial velocity")
    axs[1, 0].set_xlabel("YSOs ID")
    axs[1, 0].set_ylabel("Velocity (km/s)")
    axs[1, 0].plot(np.linspace(0, df.shape[0], 2), [df[velLSR].mean()] * 2)
    axs[1, 0].plot(np.linspace(0, df.shape[0], 2), [df[velLSR].mean()] * 2,
             label='Mean = ' + str(round(df[velLSR].mean(), 3)) + ' km/s')
    axs[1, 0].plot(np.linspace(0, df.shape[0], 2), [df[velLSR].quantile(q=0.25)] * 2,
             label='1st quartile = ' + str(round(df[velLSR].quantile(0.25), 3)) + ' km/s')
    axs[1, 0].plot(np.linspace(0, df.shape[0], 2), [df[velLSR].quantile(q=0.75)] * 2,
             label='3rd quartile = ' + str(round(df[velLSR].quantile(0.75), 3)) + ' km/s')
    axs[1, 0].legend()

    ytot = df['VTOT'].sort_values(ascending=True)
    err = df['eVTOT']
    axs[1, 1].errorbar(x, ytot, yerr=err, fmt='.', color='black', ecolor='lightgray', elinewidth=3, capsize=0, markersize=2)
    axs[1, 1].set_title("Total velocity")
    axs[1, 1].set_xlabel("YSOs ID")
    axs[1, 1].set_ylabel("Velocity (km/s)")
    axs[1, 1].plot(np.linspace(0, df.shape[0], 2), [df['VTOT'].mean()] * 2)
    axs[1, 1].plot(np.linspace(0, df.shape[0], 2), [df['VTOT'].mean()] * 2,
             label='Mean = ' + str(round(df['VTOT'].mean(), 3)) + ' km/s')
    axs[1, 1].plot(np.linspace(0, df.shape[0], 2), [df['VTOT'].quantile(q=0.25)] * 2,
             label='1st quartile = ' + str(round(df['VTOT'].quantile(0.25), 3)) + ' km/s')
    axs[1, 1].plot(np.linspace(0, df.shape[0], 2), [df['VTOT'].quantile(q=0.75)] * 2,
             label='3rd quartile = ' + str(round(df['VTOT'].quantile(0.75), 3)) + ' km/s')
    axs[1, 1].legend()

    fig.suptitle("YSO's velocities")
    plt.show()
