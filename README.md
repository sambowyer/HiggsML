# HiggsML
## An R package to accompany the SC1 group project.

This package provides some useful functions for analysing the data given in the HiggsML challenge organised by CERN. In particular, functions exist to calculate the approximate median significance (AMS) given lists of predicted and true classes (both with and without normalisation via sample weighting) and to calculate various extra features not included in the original dataset. Since our analysis of the challenge's dataset involved constructing subsets of the dataset with only positive and unlabelled (PU) data, this package also provides the functionality to generate the necessary PU labels.

NOTE: We use the `roxytest` paclage to perform testing on this packge. If you want to add tests in the roxygen details above a function via the `@tests` tag, install the `roxytest` package first.
