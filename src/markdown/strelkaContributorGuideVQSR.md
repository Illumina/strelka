<link rel='stylesheet' href='userGuide.css' />

Strelka Contributor Guide: VQSR training
========================================

Version: @WORKFLOW_VERSION@

<script src="tableOfContents.js"></script>

## Introduction

This document outlines the process of training a VQSR model for Strelka.

## Requirements

Strelka VQSR training has additional dependencies which are not included
in the starka build system.

You must have the following Python packages installed:

* numpy
* pandas
* scikit-learn
* matplotlib

## Step 1: Evaluating against a truth list

Given a VCF and a truth list, we can produce a CSV file with features and
TP/FP/FN annotation like this:

```
python ${STARKA_INSTALL}/vqsr/bin/vqsr_vcf_to_csv.py \
    src/demo/data/strelka_admix_snvs.vcf.gz \
    -o admix_training_data.csv \
    --feature-table=strelka.snv \
    --truth src/demo/data/PG_admix_truth_snvs.vcf.gz
```

## Step 2: Training a VQSR model

The next command line learns a model given a training dataset. The parameter
value of `-f` should be the same as the one used above in `--feature-table`.

```
python ${STARKA_INSTALL}/bin/vqsr_learn.py \
    -m strelka.rf \
    -f strelka.snv \
    -o model.pickle \
    admix_training_data.csv \
    --plots
```

=>

```
Reading /home/peter/workspace_lx/starka_build/admix_training_data.csv
Using default parameters.
Feature ranking:
1. feature 0:NT_REF (0.321115 +- 0.299492)
2. feature 9:N_TIER1_ALT_RATE (0.288852 +- 0.307693)
3. feature 8:T_TIER1_ALT_RATE (0.114278 +- 0.107255)
4. feature 1:QSS_NT (0.085156 +- 0.110212)
5. feature 11:MQ_ZERO_RATE (0.035398 +- 0.043017)
6. feature 7:T_DP_RATE (0.034063 +- 0.039814)
7. feature 14:SNVSB (0.033583 +- 0.074959)
8. feature 6:N_DP_RATE (0.027155 +- 0.036797)
9. feature 15:ReadPosRankSum (0.025266 +- 0.041639)
10. feature 10:MQ_SCORE (0.019818 +- 0.027179)
11. feature 3:T_FDP_RATE (0.007780 +- 0.003133)
12. feature 2:N_FDP_RATE (0.005160 +- 0.003502)
13. feature 5:T_SDP_RATE (0.001232 +- 0.000659)
14. feature 4:N_SDP_RATE (0.001144 +- 0.000810)
15. feature 13:PNOISE2 (0.000000 +- 0.000000)
16. feature 12:PNOISE (0.000000 +- 0.000000)
```

## Step 3: Calculate Quality Scores

Here, we take a given set of TPs / FPs (we use the training set in this example,
but in a real-world scenario an independent strelka run / different subsample
should be used), and write VQSR scores given the model we trained in the
previous step.

```
python ${STARKA_INSTALL}/bin/vqsr_evaluate.py -m strelka.rf -f strelka.snv \
    -c model.pickle \
    -o admix_classified.csv \
    admix_training_data.csv
```

=>

```
Reading /home/peter/workspace_lx/starka_build/admix_training_data.csv
ptag   FN     FP     TP
tag
FN    662      0      0
FP      0  32189     17
TP      0      1  15572

[3 rows x 3 columns]
```

## Step 4: Evaluate Precision / Recall for the model

```
python ${STARKA_INSTALL}/bin/vqsr_pr.py \
     -q QSS_NT,qual \
     -o admix_precisionrecall.csv \
     admix_classified.csv
```

=>

```
Reading /home/peter/workspace_lx/starka_build/admix_classified.csv
Processed 10 / 319 qual values for QSS_NT
Processed 20 / 319 qual values for QSS_NT
Processed 30 / 319 qual values for QSS_NT
Processed 40 / 319 qual values for QSS_NT
Processed 50 / 319 qual values for QSS_NT
Processed 60 / 319 qual values for QSS_NT
Processed 70 / 319 qual values for QSS_NT
Processed 80 / 319 qual values for QSS_NT
Processed 90 / 319 qual values for QSS_NT
Processed 100 / 319 qual values for QSS_NT
...
```

We can look at the result e.g. using R:

```R
data = read.csv('~/workspace_lx/starka_build/admix_precisionrecall.csv')
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
ggsave("vqsr_test.png", width=4, height=3, dpi=120)
```

![vqsr_test.png](vqsr_test.png)

(clearly, we have overtrained a bit here).

## Step 5: Export the Model for use in Strelka

Strelka uses models in JSON format:

```
python ${STARKA_INSTALL}/bin/vqsr_exportmodel.py \
    -m strelka.rf -c model.pickle -o model.json
```
