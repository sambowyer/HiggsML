# HiggsML
## An R package to accompany the SC1 group project

This package provides some useful functions for analysing the data given in the HiggsML challenge organised by CERN. In particular, functions exist to calculate the approximate median significance (AMS) given lists of predicted and true classes (both with and without normalisation via sample weighting) and to calculate various extra features not included in the original dataset. Since our analysis of the challenge's dataset involved constructing subsets of the dataset with only positive and unlabelled (PU) data, this package also provides the functionality to generate the necessary PU labels.

NOTE: We use the `roxytest` package to perform testing on this package. If you want to add tests in the roxygen details above a function via the `@tests` tag, install the `roxytest` package first.

## The Files
The package contains three R files:

- [add_pu_labels.R](R/add_pu_labels.R): In our analysis of the dataset, we create various versions of the dataset with some of the positive (background) samples mislabelled, turning the task into one of positive and unlabelled (PU) classification. This function adds columns to a given dataframe containing these partially-incorrect PU labels.
- [create_new_features.R](R/create_new_features.R): This contains only one function for ease of use which will add 70 new features to the dataset, each of which having some potentially useful physical meaning. In particular we use the formulae given by [Tim Salimans' HiggsML Challenge Github project](https://github.com/TimSalimans/HiggsML/) which came second in the original Kaggle challenge.
- [calculate_ams.R](R/calculate_ams.R): This contains functions to calculate the approximate median significance of a prediction of class labels (given the true labels) either with or without the challenge's class weights.
Weights are used to scale each data point such that the whole dataset corresponds to LHC 2012 running. Hence, if a subset $S'$ is defined, for example for testing, the weights should be renormalized:
$$w_j' = w_j \frac{\sum_i w_i \mathbb{1}\{y_i =y_j\}}{\sum_{i \in s'}w_i \mathbb{1}\{y_i =y_j\}}$$
where $y_i$ is the label of event i. 
The AMS is calculated (as per the Challenge document) via:
$$AMS = \sqrt{2\left((s+b+b_\text{reg}) \ln \left(1 + \frac{s}{b+b_\text{reg}}\right) -s \right)}$$
where $s$ and $b$ are the (potentially weighted) counts of signal and background samples respectively, and $b_\text{reg}$ is a regularization term, with a suggested value of 10 given in the Challenge.
