from astropy.io import fits
from math import *
from Velocity import *
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

#=======================================================================================================================

# convertit les coordonnées équatoriales d'ascension droite en pixel
## cube : cube de vitesses radiales
## ra : coordonnée équatoriale d'ascension droite en degrés
### rapix : coordonnée équatoriale d'ascension droite en pixel

def convertX(cube, ra):
    x0 = cube.header['CRPIX1']  # pixel de ref en X
    deltax = cube.header['CDELT1'] * 3600.  # taille des pixels en X
    rapix = x0 + ra / deltax
    return rapix

#=======================================================================================================================

# convertit les coordonnéess équatoriales de déclinaison en pixel
## cube : cube de vitesses radiales
## dec : coordonnée équatoriale de déclinaison en degrés
### decpix : coordonnée équatoriale de déclinaison en pixel

def convertY(cube, dec):
    y0 = cube.header['CRPIX2']  # pixel de ref en Y
    deltay = cube.header['CDELT2'] * 3600.  # taille des pixels en Y
    decpix = y0 + dec / deltay
    return decpix

#=======================================================================================================================

# Convertit des coordonées équatoriales RA-DEC en coordonées pixels X-Y
## cube : cube de vitesses radiales
## radec : tuple de coordonnées (ra, dec) en degrés
### coords : tuple de coordonées (ra, dec) en pixels

def convertXY(cube, radec):
    coords = (convertX(cube, radec[0]), convertY(cube, radec[1]))
    return coords

#=======================================================================================================================

# Convertit des coordonées équatoriales RA-DEC en coordonées pixels X-Y avec rotation
## cube : cube de vitesses radiales
## RA : coordonnée équatoriale d'ascension droite en degrés
## DEC : coordonnée équatoriale de déclinaison en degrés
### pix : tuple de coordonées (ra, dec) en pixels

def convertXY_abs(cube, RA, DEC):
    dec0 = cube.header['CRVAL2']  # coordonnees de reference en Y
    ra0 = cube.header['CRVAL1']  # coordonnees de reference en X
    theta = 0
    if 'CROTA1' in cube.header.keys():  # angle
        theta = cube.header['CROTA1'] * np.pi / 180.
    relra = RA - ra0  # degrees
    reldec = DEC - dec0  # degrees
    ra = (np.cos(theta) * relra + np.sin(theta) * reldec) * 3600  # arcsec
    dec = (-np.sin(theta) * relra + np.cos(theta) * reldec) * 3600  # arcsec
    pix = convertXY(cube, (ra, dec))
    return pix

#=======================================================================================================================

# Affiche une carte choisie et les YSOs du catalogue de travail selon le paramètre choisi
## df : Dataframe de travail
## plottype : 'scatter' pour nuage de points ou 'quiver' pour nuage de vecteurs
## maptype : 'Integrated Intensity' ou 'Mean Velocity' ou 'Velocity Dispersion' ou 'Lombardi'
## region : 'Main' ou 'B9' ou 'Cloak' pour les cartes de Gaudel ou 'Orion' pour la carte de Nishimura
## p : paramètre des YSOs à afficher, 'Distance' ou 'Radial Velocity' ou 'Extinction'
## velocity : nom du paramètre vitesse radiale
## zoom : 0 sans zoom, 1 avec zoom
## axis : si avec zoom, liste délimitante des axes [xmin, xmax, ymin, ymax] en pixels
## vmin : borne minimale pour le code couleur d'affichage du paramètre p
## vmax : borne maximale pour le code couleur d'affichage du paramètre p

def mapsplot(df, plottype, maptype, region, p, velocity, zoom, axis, vmin, vmax):
    # création de la figure
    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    ax.yaxis.tick_left()
    ax.yaxis.set_label_position('left')
    ax.xaxis.tick_bottom()
    ax.xaxis.set_label_position('bottom')
    ax.set_xlabel('$\delta x$ ' + "(pix)", fontsize=14)
    ax.set_ylabel('$\delta y$ ' + "(pix)", fontsize=14)
    ax.tick_params(labelsize=14)
    if zoom == 1:
        plt.xlim(axis[0], axis[1])
        plt.ylim(axis[2], axis[3])

    # filtrage du catalogue
    if p == 'Distance':
        df = blank_cleaning(df, ['ProbDist'])
        v = df['ProbDist']
        cmap = 'summer'
        title = 'Probable Distance '
        cbarLegend = 'Probable Distance of YSOs (pc)'

    if p == 'Radial Velocity':
        df = vLSR(df, velocity)
        velLSR = velocity+'LSR'
        df = blank_cleaning(df, [velLSR])
        v = df[velLSR]
        cmap = 'coolwarm'
        title = 'Radial Velocities '
        cbarLegend = 'Radial Velocities of YSOs (km/s)'

    if p == 'Extinction':
        df = blank_cleaning(df, ['Av'])
        df['NH'] = 1.8 * 10 ** 21 * df['Av']
        v = df['NH']
        cmap = 'copper'
        title = 'NH2 Quantity along sight '
        cbarLegend = 'NH2 Quantity in front of YSOs (NH/cm-2)'


    if plottype == 'quiver':
        df = blank_cleaning(df, ['pmra, pmdec'])
        title = title + 'and proper motions '


    # choix de la carte
    if maptype == 'Integrated Intensity':
        if region == 'Main':
            hdu = fits.open('intensity-map-13co10-s1-fts-Tpeak-3sigma-layer-ngc2023-24.fits')
            vminm = 0
            vmaxm = 20
        if region == 'B9':
            hdu = fits.open('intensity-map-13co10-s1-fts-Tpeak-3sigma-layer-orion-b9.fits')
            vminm = 0
            vmaxm = 12
        if region == 'Cloak':
            hdu = fits.open('intensity-map-13co10-s1-fts-Tpeak-3sigma-layer-cloak.fits')
            vminm = 0
            vmaxm = 12
        cube = hdu[0]
        data = hdu[0].data
        data = data.astype('float')
        data[data == 0.] = 'nan'
        cmapm = 'plasma'
        ec = 'black'
        cbarLegendm = 'NH2 Quantity of the cloud (NH2/cm2)'

    if maptype == 'Mean Velocity':
        if region == 'Main':
            hdu = fits.open('velocity-map-13co10-s1-fts-Tpeak-3sigma-layer-ngc2023-24.fits')
            vminm = 7
            vmaxm = 13
        if region == 'Cloak':
            hdu = fits.open('velocity-map-13co10-s1-fts-Tpeak-3sigma-layer-cloak.fits')
            vminm = 4
            vmaxm = 8
        if region == 'B9':
            hdu = fits.open('velocity-map-13co10-s1-fts-Tpeak-3sigma-layer-orion-b9.fits')
            vminm = 0
            vmaxm = 4
        if region == 'Orion':
            hdu = fits.open('Orion.CO1321.Osaka.beam204.mom1.fits')
            vminm = 1
            vmaxm = 16
        cube = hdu[0]
        data = cube.data
        cmapm = 'coolwarm'
        ec = 'black'
        cbarLegendm = 'Mean Radial Velocity of the Cloud (km/s)'

    if maptype == 'Velocity Dispersion':
        if region == 'Main':
            hdu = fits.open('dispersion-map-13co10-s1-fts-Tpeak-3sigma-layer-ngc2023-24.fits')
            vminm = 0.5
            vmaxm = 4
        if region == 'B9':
            hdu = fits.open('intensity-map-13co10-s1-fts-Tpeak-3sigma-layer-orion-b9.fits')
            vminm = 0.5
            vmaxm = 4
        if region == 'Cloak':
            hdu = fits.open('intensity-map-13co10-s1-fts-Tpeak-3sigma-layer-cloak.fits')
            vminm = 0.5
            vmaxm = 4
        cube = hdu[0]
        data = cube.data
        data = data.astype('float')
        data[data == 0.] = 'nan'
        cmapm = 'cubehelix'
        ec = 'white'
        cbarLegendm = 'Radial Velocity Dispersion of the cloud (km/s)'

    if maptype == 'Lombardi':
        hdu = fits.open('lombardi-column-density-feb18.fits')
        cube = hdu[0]
        data = cube.data
        data = data.astype('float')
        data[data == 0.] = 'nan'
        cmapm = 'copper'
        ec = 'black'
        cbarLegendm = 'Lombardi column density map'
        vminm = 10**20
        vmaxm = 10**22

    # Gestion des titres
    if plottype == 'Quiver':
        title = title + 'and proper motions '
    title = title + 'of YSOS '

    if region == 'Main':
        title = title + 'in Main Region'
    if region == 'Cloak':
        title = title + 'in Cloak Region'
    if region == 'B9':
        title = title + 'in Orion-B9 Region'
    if region == 'Orion':
        title = title + 'in Orion'
    plt.title(title, y=1.04)

    # Affichage de la carte
    im = ax.imshow(data, cmap=cmapm, vmin=vminm, vmax=vmaxm, origin='lower')  # pour plotter la carte

    # Affichage des objets
    x = df['RA']
    y = df['DEC']
    if region == 'Orion':
        df = Gal_implement(df, 'RA', 'DEC')
        x = df['glon']
        y = df['glat']
    (X, Y) = convertXY_abs(cube, x, y)

    if plottype == 'Scatter':
        cm = plt.get_cmap(cmap)
        cNorm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
        axx = plt.gca()
        divider = make_axes_locatable(axx)
        cax = divider.append_axes('right', size="5%", pad=0.1)
        caxm = divider.append_axes('right', size="5%", pad=0.8)
        cbar = plt.colorbar(scalarMap, cax=cax)
        cbarm = plt.colorbar(im, cax=caxm)
        cbar.set_label(cbarLegend, rotation=270, labelpad=15)
        cbarm.set_label(cbarLegendm, rotation=270, labelpad=15)
        ax.scatter(X, Y, c=scalarMap.to_rgba(v), vmin=vmin, vmax=vmax, marker='D', edgecolor=ec, linewidth=0.3, s=10)

    if plottype == 'Quiver':
        U, V = -df['pmra'], df['pmdec']
        if region != 'Orion':
            if 'CROTA1' in cube.header.keys():  # angle
                theta = cube.header['CROTA1'] * np.pi / 180.
                pmra = (np.cos(theta) * U + np.sin(theta) * V)  # arcsec
                pmdec = (-np.sin(theta) * U + np.cos(theta) * V)  # arcsec
                U, V = pmra, pmdec

        colors = v
        cm = plt.get_cmap(cmap)
        cNorm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
        axx = plt.gca()
        divider = make_axes_locatable(axx)
        cax = divider.append_axes('right', size="5%", pad=0.1)
        caxm = divider.append_axes('right', size="5%", pad=0.8)
        cbar = plt.colorbar(scalarMap, cax=cax)
        cbarm = plt.colorbar(im, cax=caxm)
        cbar.set_label(cbarLegend, rotation=270, labelpad=15)
        cbarm.set_label(cbarLegendm, rotation=270, labelpad=15)

        ax.quiver(X, Y, U, V, color=cm(cNorm(colors)), units='xy', scale=0.15, edgecolor=ec, linewidth=0.3)

    plt.show()

#=======================================================================================================================

# Affiche les cartes et graphiques de différences de vitesses radiales
## df : Dataframe de travail
## velocity : nom du paramètre vitesse radiale
## maptype : 'Graph' ou 'MapDiff' ou 'Correlation'
## region : 'Main' ou 'B9' ou 'Cloak' pour les cartes de Gaudel ou 'Orion' pour la carte de Nishimura
## plottype : 'scatter' pour nuage de points ou 'quiver' pour nuage de vecteurs
## type : 'abs' pour les valeurs absolues, 0 sinon.
## vmin : borne minimale pour le code couleur d'affichage du paramètre p
## vmax : borne maximale pour le code couleur d'affichage du paramètre p
## mapplot : 'Lombardi' pour afficher les YSOs sur la carte de Lombardi et 0 sinon
## thresh : Seuil de différences à partir duquel on souhaite afficher les points sur la carte

def comparison(df, velocity, maptype, region, plottype, type, vmin, vmax, mapplot, thresh):
    if maptype == 'Lombardi':
        df = blank_cleaning(df, ['Av'])
        df['NH'] = 1.8 * 10 ** 21 * df['Av']

        cmapm = 'copper'
        ec = 'grey'
        cbarLegend = 'Lombardi column density map'

        hdu = fits.open('lombardi-column-density-feb18.fits')
        vminm = 10**20
        vmaxm = 10**22

    if maptype == 'Mean Velocity':
        df = vLSR(df, velocity)
        velLSR = velocity+'LSR'
        df = blank_cleaning(df, [velLSR])

        cmapm = 'coolwarm'
        ec = 'black'
        cbarLegend = 'Radial Velocity Difference (km/s)'

        if region == 'Main':
            hdu = fits.open('velocity-map-13co10-s1-fts-Tpeak-3sigma-layer-ngc2023-24.fits')
            vminm = 7
            vmaxm = 13
        if region == 'Cloak':
            hdu = fits.open('velocity-map-13co10-s1-fts-Tpeak-3sigma-layer-cloak.fits')
            vminm = 4
            vmaxm = 8
        if region == 'B9':
            hdu = fits.open('velocity-map-13co10-s1-fts-Tpeak-3sigma-layer-orion-b9.fits')
            vminm = 0
            vmaxm = 4
        if region == 'Orion':
            hdu = fits.open('Orion.CO1321.Osaka.beam204.mom1.fits')
            vminm = 1
            vmaxm = 16

    cube = hdu[0]
    data = hdu[0].data

    df.reset_index(drop=True, inplace=True)
    x = df['RA']
    y = df['DEC']

    if region == 'Orion':
        df = Gal_implement(df, 'RA', 'DEC')
        x = df['glon']
        y = df['glat']

    (X, Y) = convertXY_abs(cube, x, y)
    A, B, V = [], [], []

    for i in range(len(X)):
        A.append(int(floor(X[i])))
        B.append(int(floor(Y[i])))
        V.append(data[B[i], A[i]])

    df['RApix'] = A
    df['DECpix'] = B
    df['Vpix'] = V
    df = df.dropna(subset=['Vpix'])
    if maptype == 'Mean Velocity':
        df['Vcomp'] = df[velLSR] - df['Vpix']
    if maptype == 'Lombardi':
        df['Vcomp'] = df['NH'] - df['Vpix']

    if type == 'abs':
        df['Vcomp'] = df['Vcomp'].abs()
        cmap = 'Greens'

    if thresh > 0:
        df = df[df['Vcomp'] > thresh]

    if thresh < 0:
        df = df[df['Vcomp'] < thresh]

    x = df['RA']
    y = df['DEC']

    if region == 'Main':
        title = 'Radial Velocity Differences between YSOs and Cloud in Main Region'
    if region == 'Cloak':
        title = 'Radial Velocity Differences between YSOs and Cloud in Cloak Region'
    if region == 'B9':
        title = 'Radial Velocity Differences between YSOs and Cloud in Orion-B9 Region'
    if region == 'Orion':
        title = 'Radial Velocity Differences between YSOs and Cloud in Orion'

    if region == 'Orion':
        df = Gal_implement(df, 'RA', 'DEC')
        x = df['glon']
        y = df['glat']
    (X, Y) = convertXY_abs(cube, x, y)

    if plottype == 'Graph':
        Vcomp = df['Vcomp'].sort_values()
        x = np.linspace(1, len(Vcomp), len(Vcomp))
        plt.xlabel('YSOs ID')
        plt.ylabel('Radial Velocity Difference (km/s)')
        plt.title(title)
        plt.plot(np.linspace(0, df.shape[0], 2), [df['Vcomp'].mean()] * 2, label='Mean = '+str(round(df['Vcomp'].mean(), 3)) + ' km/s')
        plt.plot(np.linspace(0, df.shape[0], 2), [df['Vcomp'].quantile(0.25)] * 2, label='1st quartile = '+str(round(df['Vcomp'].quantile(0.25), 3)) + ' km/s')
        plt.plot(np.linspace(0, df.shape[0], 2), [df['Vcomp'].quantile(0.75)] * 2, label='3rd quartile = '+str(round(df['Vcomp'].quantile(0.75), 3)) + ' km/s')
        plt.scatter(x, Vcomp, c='black', s=2)
        plt.legend()
        plt.show()

    if plottype == 'MapDiff':
        if mapplot == 'Lombardi':
            hdu = fits.open('lombardi-column-density-feb18.fits')
            data = hdu[0].data
            cmapm = 'copper'
            vminm = 10**20
            vmaxm = 10**22
            cbarLegendm = 'Lombardi column density map'
        cbarLegendm = 'Cloud Mean Radial Velocity (km/s) '
        fig = plt.figure()
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        ax.yaxis.tick_left()
        ax.yaxis.set_label_position('left')
        ax.xaxis.tick_bottom()
        ax.xaxis.set_label_position('bottom')
        ax.set_xlabel('$\delta x$ ' + "(pix)", fontsize=14)
        ax.set_ylabel('$\delta y$ ' + "(pix)", fontsize=14)

        im = ax.imshow(data, cmap=cmapm, vmin=vminm, vmax=vmaxm, origin='lower')  # pour plotter la carte
        v = df['Vcomp']
        if type != 'abs':
            cmap = 'PRGn'

        cm = plt.get_cmap(cmap)
        cNorm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
        if region != 'Orion':
            plt.axis([700, 1000, 100, 650])
        plt.title(title)
        ax.scatter(X, Y, c=scalarMap.to_rgba(v), vmin=vmin, vmax=vmax, marker='D', edgecolor=ec, linewidth=0.3, s=10)
        axx = plt.gca()
        divider = make_axes_locatable(axx)
        cax = divider.append_axes('right', size="5%", pad=0.1)
        caxm = divider.append_axes('right', size="5%", pad=0.8)
        cbar = plt.colorbar(scalarMap, cax=cax)
        cbarm = plt.colorbar(im, cax=caxm)
        cbar.set_label(cbarLegend, rotation=270, labelpad=15)
        cbarm.set_label(cbarLegendm, rotation=270, labelpad=15)
        plt.show()

    if plottype == 'Correlation':
        plt.scatter(df['Vpix'], df['Vcomp'], s=5)
        courbe = line_fitting(df['Vpix'], df['Vcomp'], 1)
        abs = np.linspace(min(df['Vpix']), max(df['Vpix']), 100)
        plt.plot(abs, courbe[0] + courbe[1] * abs, c='orange')
        plt.title('Correlation between RV Difference and Cloud RV')
        plt.xlabel('Cloud Radial Velocity (km/s)')
        plt.ylabel('Radial Velocity Difference (km/s)')
        plt.text(2, 8, 'Correlation = 0.1202', fontsize=10, c='orange')
        print(df[['Vpix', 'Vcomp']].corr())
        plt.show()