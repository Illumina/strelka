# Somatic Empirical Variant Score (EVS) Model Training

This is a special topic of the [Strelka User Guide](strelkaUserGuide.md).


## Introduction

This document outlines the process of training a somatic Empirical Variant Score (EVS) model for Strelka.

## Requirements

Strelka somatic EVS training has additional dependencies which are not included
in the primary build system.

You must have the following Python packages installed:

* numpy
* pandas
* scikit-learn
* matplotlib

## Step 1: Evaluating against a truth list

Given a VCF and a truth list, we can produce a CSV file with features and
TP/FP/FN annotation like this:

```
python ${STRELKA_INSTALL}/share/scoringModelTraining/somatic/bin/vcf_to_feature_csv.py \
    src/demo/data/strelka_admix_snvs_chr21:1-25000000.vcf.gz \
    -o admix_training_data.csv \
    --feature-table=strelka.snv \
    --truth src/demo/data/PG_admix_truth_snvs_chr21:1-25000000.vcf.gz
```

## Step 2: Training an EVS model

The next command line learns a model given a training dataset. The parameter
value of `-f` should be the same as the one used above in `--feature-table`.

```
python ${STRELKA_INSTALL}/share/scoringModelTraining/somatic/bin/evs_learn.py \
    -m strelka.rf \
    -f strelka.snv \
    -o model.pickle \
    admix_training_data.csv \
    --plots
```

=>

```
Reading admix_training_data.csv
Using default parameters.
Feature ranking:
1. feature 0:QSS_NT (0.360130 +- 0.040911)
2. feature 6:TIER1_ALT_RATE (0.228258 +- 0.118031)
3. feature 5:N_DP_RATE (0.155632 +- 0.048948)
4. feature 9:strandBias (0.069111 +- 0.073574)
5. feature 8:n_mapq0 (0.054176 +- 0.043498)
6. feature 7:MQ (0.049967 +- 0.032699)
7. feature 10:ReadPosRankSum (0.045467 +- 0.026747)
8. feature 1:N_FDP_RATE (0.019045 +- 0.012597)
9. feature 2:T_FDP_RATE (0.016458 +- 0.005400)
10. feature 3:N_SDP_RATE (0.000937 +- 0.000709)
11. feature 4:T_SDP_RATE (0.000820 +- 0.000598)
12. feature 12:altpos (0.000000 +- 0.000000)
13. feature 11:altmap (0.000000 +- 0.000000)
```

## Step 3: Calculate Scores

Here, we take a given set of TPs / FPs (we use the training set in this example,
but in a real-world scenario an independent strelka run / different subsample
should be used), and write EVS scores given the model we trained in the
previous step.

```
python ${STRELKA_INSTALL}/share/scoringModelTraining/somatic/bin/evs_evaluate.py -m strelka.rf -f strelka.snv \
    -c model.pickle \
    -o admix_classified.csv \
    admix_training_data.csv
```

=>

```
Reading admix_training_data.csv
ptag   FN     FP    TP
tag                   
FN    178      0     0
FP      0  17632     9
TP      0      2  5437
[3 rows x 3 columns]
```

## Step 4: Evaluate Precision / Recall for the model

```
python ${STRELKA_INSTALL}/share/scoringModelTraining/somatic/bin/evs_pr.py \
     -o admix_precisionrecall.csv \
     admix_classified.csv
```

=>

```
Reading admix_classified.csv
Processed 10 / 163 qual values for qual
Processed 20 / 163 qual values for qual
Processed 30 / 163 qual values for qual
Processed 40 / 163 qual values for qual
Processed 50 / 163 qual values for qual
Processed 60 / 163 qual values for qual
Processed 70 / 163 qual values for qual
Processed 80 / 163 qual values for qual
Processed 90 / 163 qual values for qual
Processed 100 / 163 qual values for qual
Processed 110 / 163 qual values for qual
Processed 120 / 163 qual values for qual
Processed 130 / 163 qual values for qual
Processed 140 / 163 qual values for qual
Processed 150 / 163 qual values for qual
Processed 160 / 163 qual values for qual
```

We can look at the result e.g. using R:

```R
data = read.csv('admix_precisionrecall.csv')
head(data)
```

=>

```
  X  field qual    tp   fp  fn tp_filtered fp_filtered precision    recall
1 0 QSS_NT    0 15487 1920 748           1       30286 0.8896995 0.9539267
2 1 QSS_NT    1 15460 1765 775          28       30441 0.8975327 0.9522636
3 2 QSS_NT    2 15440 1687 795          48       30519 0.9015006 0.9510317
4 3 QSS_NT    3 15422 1639 813          66       30567 0.9039329 0.9499230
5 4 QSS_NT    4 15410 1612 825          78       30594 0.9052990 0.9491839
6 5 QSS_NT    5 15403 1585 832          85       30621 0.9066988 0.9487527
```

... or make a plot like this:

```R
library(ggplot2)
ggplot(data, aes(x=recall, y=precision, color=field)) +
    geom_point() + theme_bw()
ggsave("evs_test.png", width=4, height=3, dpi=120)
```

## Step 5: Export the Model for use in Strelka

Strelka uses models in JSON format:

```
python ${STRELKA_INSTALL}/share/scoringModelTraining/somatic/bin/evs_exportmodel.py \
    -m strelka.rf -c model.pickle -o model.json
```
