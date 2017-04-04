########## 
#' research theme
#' 
#' @export
#' @name theme.research
#' @author Meridith Bartley
#' @title Research themes for ggplot2
#' @examples
#' theme_tech(theme='ants')

theme.research <- function(theme = "ants", tech_key = list(ants = list(family_title = "Circular Air Bold", 
    family_text = "Circular Air Medium", colour_title = "#674528", 
    colour_text = "#535353"))) {
    
    theme_classic() + theme(text = element_text(size = 18, family = tech_key[[theme]]$family_text)) + 
        theme(legend.title = element_blank()) + theme(plot.title = element_text(size = 25, 
        colour = tech_key[[theme]]$colour_title, family = tech_key[[theme]]$family_title)) + 
        # theme(plot.subtitle = element_text(size = 15, colour =
    # tech_key[[theme]]$colour_title,
    # family=tech_key[[theme]]$family_title)) +
    theme(axis.text.x = element_text(color = tech_key[[theme]]$colour_text)) + 
        theme(axis.text.y = element_text(color = tech_key[[theme]]$colour_text)) + 
        theme(axis.title.x = element_text(color = tech_key[[theme]]$colour_text, 
            vjust = 0)) + theme(axis.title.y = element_text(color = tech_key[[theme]]$colour_text, 
        vjust = 1.25)) + theme(plot.background = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
        theme(panel.border = element_blank()) + theme(axis.line.x = element_line(color = tech_key[[theme]]$colour_text, 
        size = 0.5), axis.line.y = element_line(color = tech_key[[theme]]$colour_text, 
        size = 0.5)) + theme(line = element_line(color = tech_key[[theme]]$colour_text)) + 
        theme(rect = element_rect(color = tech_key[[theme]]$colour_text)) + 
        theme(axis.ticks.x = element_line(color = tech_key[[theme]]$colour_text)) + 
        theme(axis.ticks.y = element_line(color = tech_key[[theme]]$colour_text))
}
