#' Add PU Labels
#'
#' To simulate positive and unlabelled (PU) data we denote the background samples as positive with label 1 and signal samples as negative with label 0, before misclassifying some proportion (1-pu.pi) of the background samples.
#' This function adds the resulting labels as columns in the dataset it is passed as an argument.
#'
#' @param df A dataframe containing data from the HiggsML challenge (must have true labels of "s" and "b" for signal and background stored in a column named "Label").
#' @param pu.pi The proportion of background samples to maintain their correct labels.
#'
#' @return The dataframe df with an extra column containing the PU labels. The column will be named "pu[pu.pi]" where [pu.pi] represents the numeric value of [pu.pi] as a string.
#' @export
#'
#' @examples
#' \dontrun{
#'   add_pu_labels(training_dataset, 0.6)
#' }
#'
#' @tests
#' allBackground = data.frame(Label = rep("b", 100))
#' allOnes = rep(1, 100)
#' expect_equal(add_pu_labels(allBackground, 1)$pu1, allOnes)
#'
#' allSignal= data.frame(Label = rep("b", 100))
#' allZeroes = rep(0, 100)
#' expect_equal(add_pu_labels(allSignal, 1)$pu1, allOnes)
#' expect_equal(add_pu_labels(allBackground, 0)$pu0, allZeroes)
#'
#' expect_equal(sum(add_pu_labels(allBackground, 0.5)$pu0.5), 50)
#'
#' df = read.csv("../../data/atlas-higgs-challenge-2014-v2.csv")
#' trueBackgroundCount = sum(df$Label == "b")
#' for (pu.pi in c(seq(0, 0.8, 0.2), 1)){
#'   df = add_pu_labels(df, pu.pi)
#' }
#' expect_equal(sum(df$pu0), 0)
#' expect_equal(sum(df$pu0.2), ceiling(trueBackgroundCount*0.2))
#' expect_equal(sum(df$pu0.4), ceiling(trueBackgroundCount*0.4))
#' expect_equal(sum(df$pu0.6), ceiling(trueBackgroundCount*0.6))
#' expect_equal(sum(df$pu0.8), ceiling(trueBackgroundCount*0.8))
#' expect_equal(sum(df$pu1), trueBackgroundCount)
add_pu_labels <- function(df, pu.pi){
  if (pu.pi < 1){
    pu_label <- rep(1, nrow(df))
    pu_label[df$Label == "s"] <- 0  # make all signal samples unlabelled (i.e. have label 0)

    #Mislabel (1-pu.pi) proportion of the background samples
    pu_label[sample(which(df$Label=="b"), (1-pu.pi)*sum(df$Label=="b"), replace=F)] <- 0

    # Add the pu labels to the dataframe
    df = cbind(df,pu_label)
    colnames(df)[ncol(df)] <- paste("pu", pu.pi, sep="")

  } else if (pu.pi == 1) {
    # Separate case for pu.pi==1 to avoid potential sampling bugs
    df$pu1 = as.numeric(df$Label=="b")
  }

  df
}
