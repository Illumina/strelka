# Somatic Empirical Variant Score (EVS) Model Training

This is a special topic of the [Strelka User Guide](strelkaUserGuide.md).

## Introduction

This document outlines the process of training a somatic Empirical Variant Score (EVS) model for Strelka.

## Requirements

Strelka somatic EVS training has additional dependencies which are not included
in the primary build system. All packages required to retrain the EVS model and
run the steps in this guide are provided in the [training environment
Dockerfile](trainingSomaticEmpriicalScore/Dockerfile), which can be used to
either setup an EVS training docker image or as documentation to install
all dependencies on another system.

## Step 1: Build a training data set

For EVS training, the strelka workflow must be configured with the optional
`--reportEVSFeatures` argument. This will add a new VCF INFO field called `EVSF`
to both SNV and indel VCF outputs. EVS features will be reported with this flag even
when scoring itself is turned off with the `--disableEVS` option.

### Step 1a: Evaluate a single strelka somatic analysis against a truth list

Given a VCF and a truth list, a CSV file with features and
TP/FP/FN annotation can be produced as follows:

```
python ${STRELKA_INSTALL_PATH}/share/scoringModelTraining/somatic/bin/vcf_to_feature_csv.py \
    --truth ${STRELKA_INSTALL_PATH}/share/demo/strelka/data/PG_admix_truth_snvs_chr21:1-25000000.vcf.gz \
    --features strelka.snv \
    --output admix_training_data.csv \
    ${STRELKA_INSTALL_PATH}/share/demo/strelka/data/strelka_admix_snvs_chr21:1-25000000.vcf.gz
```

### Step 1b: Handle multiple analyses

When multiple analyses are available for a training procedure, process each analysis
using the training procedure described in Step 1a, keeping a separate csv output file
for each. These training data files can be combined as required for the training
procedure described below.


## Step 2: Training an EVS model

The next step is to train a model given one or more training datasets.
The `--features` argument below must match that used to generate
the input training dataset(s).

```
python ${STRELKA_INSTALL_PATH}/share/scoringModelTraining/somatic/bin/evs_learn.py \
    --features strelka.snv \
    --model strelka.rf \
    --plots \
    --output admix_training_model.pickle \
    admix_training_data.csv
```

=>

```
Reading admix_training_data.csv
Using default parameters.
Feature ranking:
1. feature 0:QSS_NT (0.341969 +- 0.040246)
2. feature 6:TIER1_ALT_RATE (0.190468 +- 0.123727)
3. feature 5:N_DP_RATE (0.154801 +- 0.068986)
4. feature 11:altmap (0.058069 +- 0.062328)
5. feature 8:n_mapq0 (0.056272 +- 0.051838)
6. feature 7:MQ (0.051763 +- 0.043005)
7. feature 12:altpos (0.047467 +- 0.053393)
8. feature 10:ReadPosRankSum (0.032803 +- 0.020241)
9. feature 9:strandBias (0.032238 +- 0.056309)
10. feature 1:N_FDP_RATE (0.018251 +- 0.013008)
11. feature 2:T_FDP_RATE (0.014174 +- 0.007238)
12. feature 4:T_SDP_RATE (0.000968 +- 0.000762)
13. feature 3:N_SDP_RATE (0.000757 +- 0.000641)
```

## Step 3: Calculate Scores

Given a trained model any set of TP/FP variants can be scored to evaluated the trained
model. In this example the training data itself is used is used for simplicity, although
independent test data should be used for more accurate evaluation in a real training
scenario.

```
python ${STRELKA_INSTALL_PATH}/share/scoringModelTraining/somatic/bin/evs_evaluate.py \
    --features strelka.snv \
    --classifier admix_training_model.pickle \
    --output admix_classified.csv \
    admix_training_data.csv
```

=>

```
Reading admix_training_data.csv
ptag   FN     FP    TP
tag                   
FN    178      0     0
FP      0  17640     1
TP      0      0  5439

[3 rows x 3 columns]
```

## Step 4: Evaluate Precision / Recall for the model

```
python ${STRELKA_INSTALL_PATH}/share/scoringModelTraining/somatic/bin/evs_pr.py \
     --output admix_precisionrecall.csv \
     admix_classified.csv
```

=>

```
Reading admix_classified.csv
Processed 10 / 99 qual values for qual
Processed 20 / 99 qual values for qual
Processed 30 / 99 qual values for qual
Processed 40 / 99 qual values for qual
Processed 50 / 99 qual values for qual
Processed 60 / 99 qual values for qual
Processed 70 / 99 qual values for qual
Processed 80 / 99 qual values for qual
Processed 90 / 99 qual values for qual
```

We can look at the result e.g. using R:

```R
data = read.csv('admix_precisionrecall.csv')
head(data)
```

=>

```
  X field qual   tp   fp  fn tp_filtered fp_filtered precision    recall
1 0  qual 0.00 5439 5531 178           0       12110 0.4958067 0.9683105
2 1  qual 0.01 5439 3508 178           0       14133 0.6079133 0.9683105
3 2  qual 0.02 5439 2564 178           0       15077 0.6796201 0.9683105
4 3  qual 0.03 5439 1961 178           0       15680 0.7350000 0.9683105
5 4  qual 0.04 5439 1521 178           0       16120 0.7814655 0.9683105
6 5  qual 0.05 5439 1200 178           0       16441 0.8192499 0.9683105
```

... or make a plot like this:

```R
library(ggplot2)
ggplot(data, aes(x=recall, y=precision, color=field)) +
    geom_point() + theme_bw()
ggsave("admix.png", width=4, height=3, dpi=120)
```

An example plot from this procedure is:

![Training data ROC curve](trainingSomaticEmpriicalScore/admix.png)

...showing the expected result of testing on a very small training data set.

## Step 5: Export the model for use in Strelka

Strelka uses models in JSON format, which can be produced from the model pickle file as follows:

```
python ${STRELKA_INSTALL_PATH}/share/scoringModelTraining/somatic/bin/evs_exportmodel.py \
    --classifier admix_training_model.pickle \
    --output admix_training_model.json
```
