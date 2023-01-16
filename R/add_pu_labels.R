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
add_pu_labels <- function(df, pu.pi){
  if (pu.pi < 1){
    pu_label <- rep(1, nrow(df))
    pu_label[df$Label == "s"] <- 0
    pu_label[sample(which(df$Label=="b"), (1-pu.pi)*sum(df$Label=="b"), replace=F)] <- 0
    df = cbind(df,pu_label)
    colnames(df)[ncol(df)] <- paste("pu", pu.pi, sep="")
  } else if (pu.pi == 1) {
    # Separate case to avoid potential sampling bugs
    df$pu1 = as.numeric(df$Label=="b")
  }

  df
}
