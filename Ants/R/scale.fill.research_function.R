##########
#' research scale fill
#' 
#' @export
#' @name scale.fill.research
#' @author Meriidth Bartley
#' @title Research scale fill for ggplot2
#' @examples
#' scale_fill_tech(theme="ants")


scale.fill.research <- function(theme="ants", tech_key = list(
  ants = c("#bc5356", "#bcb953", "#538bbc", "#53bc84",  "#8453bc")
)) {
  
  scale_fill_manual(values=tech_key[[theme]])
  
}