# Score biologique — version k-mer

Portage du **score biologique** du dépôt [`Score-biologique`](https://github.com/ambregly/Score-biologique)
(pipeline Chromoanagenesis-AML, basé sur les sorties **Arriba**) vers une version
qui fonctionne à partir de **fichiers pauvres en information** : les fichiers **k-mer**.

Ici on **n'a plus Arriba**. Le script **reconstruit** les caractéristiques des
fusions manquantes (brin → direction → **type / classe chimérique**) à partir du
peu qu'on a, et remplace les composantes Arriba par des équivalents k-mer.

---

## Ce qui change par rapport au score original

| | Score original (Arriba) | Version k-mer |
|---|---|---|
| **type** chimérique | lu dans Arriba | 🔧 **reconstruit** (chr + brin + direction + distance) |
| **reading_frame** | prédiction peptidique Arriba | 🔧 **reconstruit** via GENCODE (bordure d'exon / CDS) |
| **spécificité** | max WT vs patho | 🔧 contraste patho ↔ blacklist normaux |
| **reads** Arriba | `split_reads1+2` | ✅ remplacé par le **comptage k-mer** (intégré à la spécificité) |
| **confidence** | high/medium/low | ❌ **supprimée** |
| **WHO** | liste WHO | ✅ conservée |
| **agrégation** | MAX par fusion | **MOYENNE** par fusion |

---

## Entrées

### 1. Fichier(s) k-mer patho — TSV, sans en-tête, 3 colonnes
```
<gene5'_chr5'_bp5'_gene3'_chr3'_bp3'_index|index...>   <patient>   <comptage>
ACTB_chr7_5570155_PTMAP5_chr13_82264248_813|2628       BA2816R     7.19
```
- Le **1er gène est en 5'**, le 2e en 3'.
- Une ligne = un k-mer ; une même fusion apparaît sur plusieurs lignes (plusieurs
  index, plusieurs patients).

### 2. Blacklist des normaux — `kmer_normaux_33.tsv`
Même format que (1), comptages k-mer chez les ~33 échantillons **normaux**. Le
nombre de normaux est lu depuis le nom du fichier (`..._33.tsv` → 33).

### 3. Deux BED — `left` (gène 5') et `right` (gène 3')
```
chr    start     end       <gene_chr_bp_index...>              .   strand
chr18  24081189  24081214  KCTD1_chr18_24081214_4|25|37|...    .   -
```
→ apportent le **brin (+/−)** de chaque gène.

### 4. Annotation GENCODE (hg38)
Téléchargée **une seule fois** puis mise en **cache** (`annot_cache/`). Fournir un
fichier local avec `--gtf` pour éviter tout accès réseau.

---

## Reconstructions

**Brin → direction** (partenaire de fusion) :

| Rôle | Brin | direction |
|---|---|---|
| 5' (BED left) | `+` | downstream |
| 5' (BED left) | `−` | upstream |
| 3' (BED right) | `+` | upstream |
| 3' (BED right) | `−` | downstream |

**Direction + chr + distance → type (classe chimérique)** :

| Situation | type | classe |
|---|---|---|
| chromosomes différents | Translocation | Classe 1 |
| même chr, directions **égales** | Inversion | Classe 4 |
| même chr, dir. opposées, délétion colinéaire même brin < seuil kb | Read-through | Classe 2 |
| même chr, dir. opposées, `bp5' < bp3'` (dir 5'=downstream) | Délétion | Classe 2 |
| même chr, dir. opposées, ordre inverse | Duplication | Classe 3 |

(Brin manquant → type `Indéterminé` ; une translocation reste détectable sans brin.)

---

## Le score

```
score_norm = Σ(fraction × poids) / Σ(poids)          (max 11 pts par défaut)
```

| Composante | Poids | Fraction ∈ [0,1] |
|---|---|---|
| **spécificité** | 4 | `⅔·(1 − présence_normaux) + ⅓·expr_patho` |
| **type** | 3 | Transloc./Inversion 1 · Délét./Duplic. 0,5 · Read-through 0,25 · sinon 0 |
| **reading_frame** | 2 | moyenne des 2 breakpoints : bordure exacte 1 · <5 nt 0,9 · CDS 0,5 · reste 0 |
| **WHO** | 2 | fusion d'intérêt WHO → 1 · sinon 0 |

> La **spécificité** a le poids le plus fort (elle prime sur le type chimérique).
> Tous les poids restent surchargeables en CLI (`--spec --type --frame --who`).

avec, sur l'ensemble des fusions candidates :
- `présence_normaux = moyenne(freq_norm, val_norm)` — l'**absence des normaux domine** ;
- `freq_norm = n_normaux_positifs / N_normaux` ;
- `val_norm  = minmax(log1p(comptage_moyen_normaux))` ;
- `expr_patho = minmax(log1p(comptage_moyen_patho))` — plus la couverture k-mer
  patho est haute, mieux c'est.

**Agrégation des comptages = MOYENNE** sur tous les k-mers × échantillons *positifs*
d'une fusion (les absents ne diluent pas). Vaut pour patho **et** normaux.

**Une fusion = une paire de gènes.** Une même paire (`gene5'--gene3'`) peut
apparaître avec plusieurs breakpoints (plusieurs k-mers) : elle n'est comptée
qu'**une seule fois**, en retenant la variante au **score le plus haut**. Toutes
les variantes restent tracées dans `fusions_score_kmer_all_variants.tsv`.

**Matching blacklist** : clé `gene5'_bp5'_gene3'_bp3'` (avec breakpoints, **sans**
chromosomes ni index).

**Priorités** : `P1 ≥ 65 %` · `P2 ≥ 40 %` · `P3 ≥ 20 %` · sinon `NP`.
Un poids à `0` retire la composante du calcul **et** du dénominateur.

---

## Cartes de priorisation & concentration

Les cartes (`carte_priorisation_*.pdf`) sont des *bubble scatter* :
**x** = score biologique, **y** = comptage k-mer **max chez un patient** (log),
**taille** = concentration, **couleur** = répartition entre patients (nb de patients
positifs) ou type chimérique. La zone ombrée = priorité P1 (score ≥ 65 %).

### Concentration (part du patient principal)

$$\text{concentration} = \frac{\max}{\text{somme}} \in\ ]0,1]$$

où `max` = comptage de la fusion **chez son patient le plus positif** et `somme` =
total de ses comptages **sur tous les patients positifs**. C'est la **part de
l'expression totale portée par le patient principal**, exprimée en **pourcentage** :

- **100 %** = toute l'expression est chez **un seul** patient → événement *privé* ;
- **petit %** = expression **étalée** sur beaucoup de patients → diffuse.

**Plus la bulle est grosse, plus la fusion est concentrée sur un seul patient** —
exactement le profil recherché pour un événement spécifique à un individu / sous-groupe.

| Comptages par patient | max | somme | concentration | Lecture |
|---|---|---|---|---|
| `[100, 0, 0]` | 100 | 100 | **100 %** | tout chez 1 patient → grosse bulle |
| `[100, 100, 100]` | 100 | 300 | **33 %** | étalée sur 3 → moyenne |
| `[10, 5, 5]` | 10 | 20 | **50 %** | moitié sur le patient principal |
| `[3, 3, 3, 3, 3]` | 3 | 15 | **20 %** | diffuse sur 5 → petite bulle |

> La concentration est **indépendante du niveau d'expression** (le niveau se lit sur
> l'axe **y**). Les deux se combinent : haut sur y **et** grosse bulle = fusion à la
> fois fortement exprimée **et** privée.

La **charge par patient** (`charge_par_patient.pdf`) affiche des **codes anonymes**
(`P01` = charge la plus forte) car les patients sont nombreux ; la correspondance
code ↔ vrai nom (avec la charge et le détail par type) est dans
`charge_par_patient_correspondance.tsv`.

## Utilisation

```bash
Rscript score_kmer.R \
    --dir-kmer  kmer_score/patho \
    --blacklist kmer_score/kmer_normaux_33.tsv \
    --bed-left  kmer_score/left.bed \
    --bed-right kmer_score/right.bed \
    --dir-out   analyse_fusions_kmer
    # --gtf annot_cache/gencode.v45.annotation.gtf.gz   # optionnel (sinon téléchargé)
```

### Options

| Option | Défaut | Rôle |
|---|---|---|
| `--dir-kmer` | `kmer_score/patho` | dossier des `*.tsv` k-mer patho (ou un fichier unique) |
| `--blacklist` | `kmer_score/kmer_normaux_33.tsv` | blacklist des normaux |
| `--bed-left` / `--bed-right` | `kmer_score/left.bed` / `right.bed` | BED brin 5' / 3' |
| `--gtf` | *(téléchargé)* | GTF GENCODE local (évite le réseau) |
| `--gtf-url` | GENCODE v45 | URL de téléchargement |
| `--cache-dir` | `annot_cache` | cache de l'annotation |
| `--dir-out` | `analyse_fusions_kmer` | dossier de sortie |
| `--spec` `--type` `--frame` `--who` | 4 / 3 / 2 / 2 | poids des composantes (0 = retirée) |
| `--n-top` | 30 | fusions affichées dans les figures |
| `--rt-kb` | 300 | seuil read-through (kb) |
| `--bp-tol` | 5 | distance max (nt) pour « proche d'une bordure d'exon » (fraction 0,9) |

### Sorties (`--dir-out`)
```
analyse_fusions_kmer/
├── fusions_score_kmer.tsv              # 1 ligne par paire de gènes (meilleur score)
├── fusions_score_kmer_all_variants.tsv # toutes les variantes de breakpoints
├── fusions_specifiques_kmer.tsv        # sous-ensemble absent des normaux
├── charge_par_patient_correspondance.tsv # code anonyme ↔ vrai nom de patient (+ charge)
└── figures/
    ├── carte_priorisation_repartition.pdf # score × expression max/patient, couleur = répartition
    ├── carte_priorisation_type.pdf     # même carte, couleur = type chimérique
    ├── charge_par_patient.pdf          # nb de fusions chromo-spé. par patient (codes anonymes), par type
    ├── score_classement.pdf            # top fusions par score
    ├── repartition_types.pdf           # types chimériques reconstruits
    └── score_decomposition.pdf         # décomposition du score par composante
```

---

## Exemple reproductible

Le dossier [`example/`](example/) contient un mini-jeu synthétique (k-mer patho,
blacklist, BED, mini-GTF) qui exerce chaque cas (translocation, inversion,
read-through, WHO, présence chez les normaux) :

```bash
Rscript score_kmer.R \
    --dir-kmer  example/kmer_score/patho \
    --blacklist example/kmer_score/kmer_normaux_33.tsv \
    --bed-left  example/kmer_score/left.bed \
    --bed-right example/kmer_score/right.bed \
    --gtf       example/mini.gtf \
    --dir-out   example/out
```

---

## Prérequis

- **R ≥ 4.1** avec `tidyverse` et `scales` (obligatoires) ; `ggrepel` (optionnel —
  labels non chevauchants sur les cartes, sinon `geom_text`).
- Outils shell `zcat`, `awk`, `grep` (pré-filtrage de l'annotation GENCODE).
- Un device **cairo** + une **locale UTF-8** (fixés automatiquement par le script)
  pour le rendu des accents et du `²` dans les figures.

## Notes méthodologiques

- Reconstruction du type inspirée de la logique Arriba (orientation des breakpoints).
- Read-through = délétion colinéaire, même brin, distance < `--rt-kb` (défaut 300 kb).
- Les normalisations `minmax`/`log1p` sont **relatives à la cohorte analysée**.
- Références : Rufflé 2017 & 2024.
