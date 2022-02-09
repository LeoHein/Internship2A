import matplotlib.pyplot as plt
import numpy as np

#=======================================================================================================================

# Affiche par ordre croissant les distances probables et erreurs associées des objets dans le Dataframe d'entrée
## df : Dataframe de travail

def distance_graph(df):
    df = df.sort_values(by='ProbDist', ascending=True)
    x = [i for i in range(len(df))]
    y = df['ProbDist']
    err = (df['MaxDist'] - df['MinDist']) / 2
    plt.errorbar(x, y, yerr=err, fmt='.', color='black', ecolor='lightgray', elinewidth=3, capsize=0, markersize = 2)

    plt.plot(np.linspace(0, df.shape[0], 2), [df['ProbDist'].mean()] * 2, label='Mean = ' + str(round(df['ProbDist'].mean(), 3)) + ' pc')
    plt.plot(np.linspace(0, df.shape[0], 2), [df['ProbDist'].quantile(q=0.25)] * 2, label='1st quartile = ' + str(round(df['ProbDist'].quantile(0.25), 3)) + ' pc')
    plt.plot(np.linspace(0, df.shape[0], 2), [df['ProbDist'].quantile(q=0.75)] * 2, label='3rd quartile = ' + str(round(df['ProbDist'].quantile(0.75), 3)) + ' pc')
    plt.xlabel("ID YSO")
    plt.ylabel("Probable Distance + error (parsec)")
    plt.title("Probable distance of the YSOs in ascending order")
    plt.legend(loc='upper left')
    plt.show()

#=======================================================================================================================

# Filtre les objets du Dataframe d'entrée selon leurs distances probables et leurs erreurs associées
## df : Dataframe de travail
## MinDist, MaxDist : respectivement bornes inférieure et supérieure de l'intervalle de distances accepté en pc
## error : 0 pour aucun filtre sur l'erreur ou un nombre positif comme seuil d'erreur maximum accepté en pc
### df : Dataframe modifié

def position_filter(df, MinDist, MaxDist, error):
    if error == 0:
        df = df[(df['ProbDist'] > MinDist) & (df['ProbDist'] < MaxDist)]
    else:
        df = df[(df['ProbDist'] > MinDist) & (df['ProbDist'] < MaxDist)]
        df['epsilon'] = (df['parallax'] + df['parallax_error']) * (df['parallax'] - df['parallax_error']) / 1000 * error
        df = df[df['parallax_error'] < df['epsilon']]
    return df
