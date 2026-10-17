# Optimisation de portefeuille sous contraintes (R)

Construction et évaluation de portefeuilles actions à partir de rendements mensuels :
minimum-variance, réplication d'indice et frontière efficiente, avec comparaison de deux
solveurs et une **évaluation hors échantillon**.

## Ce que fait le script

1. **Période d'apprentissage** (60 premiers mois), où se font toutes les estimations :
   - covariance estimée par **shrinkage de Ledoit-Wolf**. Avec plus d'actions que de mois,
     la covariance empirique est singulière et l'optimiseur exploiterait du bruit
     d'estimation ;
   - portefeuille **minimum-variance global** sous contraintes (long-only, 0,01 % ≤ wᵢ ≤ 20 %) ;
   - portefeuille de **réplication d'indice** (minimisation de la tracking error) ;
   - **frontière efficiente** résolue par deux méthodes : programmation quadratique exacte
     (`quadprog`) et optimisation non linéaire SLSQP (`nloptr`). On compare les volatilités
     obtenues et l'écart maximal de poids entre les deux.
2. **Période de test** (mois restants, jamais utilisés pour estimer quoi que ce soit) :
   - minimum-variance à poids fixes (estimés sur l'apprentissage) ;
   - minimum-variance **glissant** : chaque mois, poids ré-estimés sur tout l'historique
     disponible jusque-là ;
   - réplication d'indice, portefeuille équipondéré (1/N) et indice comme références.
3. **Métriques** (codées en R de base, sans boîte noire) : rendement et volatilité
   annualisés, Sharpe, Sortino, Omega, bêta, Treynor, alpha de Jensen, max drawdown,
   Calmar, VaR et ES historiques à 95 %, tracking error, information ratio.

## Résultats

<!-- À compléter après exécution : copier ici le tableau outputs/metrics_out_of_sample.csv
     et insérer les figures ci-dessous. -->

![Frontière efficiente](outputs/figures/efficient_frontier.png)
![Hors échantillon](outputs/figures/out_of_sample_growth.png)

## Données

Fichier Excel `data/returns.xlsx` contenant des rendements **mensuels simples en décimal** :

| Colonne | Contenu |
|---|---|
| `Date` | fin de mois |
| `Bond` | taux sans risque (unité réglable via `RF_UNIT` en tête de script) |
| `INDEX` | rendement de l'indice de référence |
| autres colonnes | rendements des actions de l'univers |

## Exécution

```r
# depuis le dossier du projet
source("portfolio_optimization.R")
```

Les paramètres (chemin des données, longueur de l'apprentissage, bornes de poids, unité
du taux sans risque) sont regroupés en tête du script. Résultats : `outputs/*.csv`,
figures : `outputs/figures/`.

## Limites

- La période de test est courte (quelques mois) : les écarts de performance entre
  stratégies ne sont pas statistiquement significatifs.
- Les rendements espérés (frontière) sont des moyennes historiques, très bruitées ;
  la frontière est illustrative. Le minimum-variance, qui n'utilise pas les rendements
  espérés, est le portefeuille le plus robuste du lot.
- Pas de coûts de transaction, alors que le portefeuille glissant change de poids chaque mois.
