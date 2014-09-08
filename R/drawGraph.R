#' Draw graph function
#' 
#' A function to draw the target inhibition network.
#' 
#' @param draw_data a data frame combines drug-target interaction data with drug sensitivity. The column names
#' must be upper case.
#' @return An image of the target inhibition network.
#' @author Jing Tang \email{jing.tang@@helsinki.fi} 
#' @references Tang J, Karhinen L, Xu T, Szwajda A, Yadav B, Wennerberg K, Aittokallio T. 
#' Target inhibition networks: predicting selective combinations of druggable targets to block cancer 
#' survival pathways. PLOS Computational Biology 2013; 9: e1003226.
#' @examples 
#' data(tyner_interaction_binary)
#' data(tyner_sensitivity)
#' x<-data.frame(tyner_interaction_binary[, 1:6])
#' y<-tyner_sensitivity[,1]
#' #binarize the sensitivity data
#' one<-which(y>0.5)
#' zero<-which(y<=0.5)
#' SENS<-y
#' SENS[one]<-1
#' SENS[zero]<-0
#' draw_data<-cbind(x, SENS)
#' drawGraph(draw_data)
drawGraph <- function(draw_data) {
    
    # column names must be upper case
    
    column_names <- dimnames(draw_data)[[2]]
    dimnames(draw_data)[[2]][length(column_names)] <- "SENS"
    boolean <- eqmcc(draw_data, outcome = "SENS")
    boolean$essential
    num_components <- length(boolean$essential)
    expressions <- sapply(boolean$essential, function(x) strsplit(x, "\\*"))
    expressions_compact <- expressions
    
    # only treat the true conditions
    is.upper <- "[A-Z]"
    for (i in 1:num_components) {
        result <- grepl(pattern = is.upper, x = expressions[[i]])
        expressions_true <- expressions[[i]][which(result == TRUE)]
       
        #expressions_compact[[i]] <- expressions[[i]][which(result == T)]
        expressions_compact[[i]] <- unlist(lapply(expressions[[i]][which(result == TRUE)],
                                                  function(x) unlist(strsplit(x, ".", fixed=TRUE))[1]))
    }
    height_components <- lapply(expressions_compact, length)
    
    seg <- 1
    seg_in <- 0.2
    margin <- 0.2
    height_figure <- max(unlist(height_components))
    width_figure <- num_components * (1 + margin)
    
    # drawing
    dummy <- 0
    plot(dummy, dummy, type = "n", axes = FALSE, ann = FALSE, xlim = c(1, width_figure + 1), ylim = c(0, height_figure))  # no axes, no labels
    lines(c(1, 1 + margin), c(height_figure/2, height_figure/2), type = "l")  # starting
    
    for (i in 1:num_components) {
        leg <- height_components[[i]]
        # y = vector(mode='numeric',length = leg)
        y0 <- height_figure/2 - (leg - 1)/2
        y0 <- seq(y0, y0 + leg - 1)
        x0 <- rep(i * (seg + margin), leg)
        x1 <- rep(i * (seg + margin) + seg_in, leg)
        y1 <- y0
        x2 <- x0 + seg
        x3 <- x0 + seg - margin
        segments(x0, y0, x1, y1)
        segments(x2, y0, x3, y1)
        lines(x0[c(1, 1)], y0[c(1, leg)], type = "l")
        lines(x2[c(1, 1)], y0[c(1, leg)], type = "l")
        lines(c(x2[1], x2[1] + margin), c(height_figure/2, height_figure/2), type = "l")
        
        # write names
        text(x1, y0, x3, y1, labels = expressions_compact[[i]], pos = 4)
    }
    
} 
