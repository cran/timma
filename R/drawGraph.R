#' Draw graph function
#' 
#' A function to draw the target inhibition network.
#' 
#' @param draw_data a data frame combines drug-target interaction data with drug sensitivity. The column names
#' must be upper case.
#' @return An image of the target inhibition network.
#' @author Jing Tang \email{jing.tang@@helsinki.fi}, Liye He \email{liye.he@@helsinki.fi}
#' @references Tang J, Karhinen L, Xu T, Szwajda A, Yadav B, Wennerberg K, Aittokallio T. 
#' Target inhibition networks: predicting selective combinations of druggable targets to block cancer 
#' survival pathways. PLOS Computational Biology 2013; 9: e1003226.
#' @examples 
#' \dontrun{
#' data(tyner_interaction_binary)
#' data(tyner_sensitivity)
#' y<-tyner_sensitivity[,1]
#' k_selected<-sffs(tyner_interaction_binary, y)$k_sel
#' x<-data.frame(tyner_interaction_binary[, k_selected])
#' #binarize the sensitivity data
#' one<-which(y>0.5)
#' zero<-which(y<=0.5)
#' SENS<-y
#' SENS[one]<-1
#' SENS[zero]<-0
#' draw_data<-cbind(x, SENS)
#' drawGraph(draw_data)
#' }

drawGraph <- function(draw_data) {
    
    # column names must be upper case
    
    column_names <- dimnames(draw_data)[[2]]
    dimnames(draw_data)[[2]][length(column_names)] <- "SENS"
    boolean <- eqmcc(draw_data, outcome = "SENS")
    boolean$essential
    num_components <- length(boolean$essential)
    expressions <- sapply(boolean$essential, function(x) strsplit(x, "\\*"))
    
    expressions_compact <- list()
    # only treat the true conditions
    is.upper <- "[A-Z]"
    num <- 1
    for (i in 1:num_components) {
        
        result <- grepl(pattern = is.upper, x = expressions[[i]])
        if(any(result)){
          expressions_compact[[num]] <- unlist(lapply(expressions[[i]][which(result == TRUE)],
                                                    function(x) unlist(strsplit(x, ".", fixed=TRUE))[1]))
          num <- num + 1
        }
        
    }
    num_components <- num - 1
    height_components <- lapply(expressions_compact, length)
    
    seg <- 1
    seg_in <- 0.2
    margin <- 0.2
    height_figure <- max(unlist(height_components))
    width_figure <- num_components * (1 + max(nchar(unlist(height_components)))*0.4)
    #width_figure <- num_components * (1 + margin)
    
    # drawing
    dummy <- 0
    pdf(file="targetInhibitionNetwork.pdf", width=12, height=12)
    plot(dummy, dummy, type = "n", axes = FALSE, ann = FALSE, xlim = c(1, width_figure + 1), ylim = c(0, height_figure))  # no axes, no labels
    lines(c(1, 1 + margin), c(height_figure/2, height_figure/2), type = "l")  # starting
    start_line <- 1 + margin
    
    for (i in 1:num_components) {
        leg <- height_components[[i]]
        len_target <- c()
        for(j in 1:leg){
          len_target <- c(len_target, nchar(expressions_compact[[i]][j]))
        }
        max_len <- max(len_target)
        # y = vector(mode='numeric',length = leg)
        y0 <- height_figure/2 - (leg - 1)/2
        y0 <- seq(y0, y0 + leg - 1)
        #x0 <- rep(i * (seg + margin), leg)
        x0 <- rep(start_line, leg)
        #x1 <- rep(i * (seg + margin) + seg_in, leg)
        x1 <- rep(x0+seg_in, leg)
        y1 <- y0
        x2 <- x1 + seg_in*max_len/4*num/3
        #x2 <- rep(1.6, leg)
        x3 <- x2 + seg_in
        #x3 <- x0 + seg - margin
        segments(x0, y0, x1, y1)
        segments(x2, y0, x3, y1)
        lines(x0[c(1, 1)], y0[c(1, leg)], type = "l")
        lines(x3[c(1, 1)], y0[c(1, leg)], type = "l")
        lines(c(x3[1], x3[1] + margin), c(height_figure/2, height_figure/2), type = "l")
        start_line <- x3[1] + margin
        
        # write names
        text(x1, y0, labels = expressions_compact[[i]], pos = 4)
    }
    
    #lines(c(start_line, start_line + margin), c(height_figure/2, height_figure/2), type = "l")  # ending
    dev.off()
    
    # two terminal version
    cat(c(),file="targetInhibitionNetwork.sif")
    Terminals <- paste("T",1:2, sep="") # the number of joints for the SIF format
    
    # connecting the first compoonent with T1
    for (j in expressions_compact[[1]]){
      lines <- paste(c(Terminals[1],"xx",j), collapse="\t")
      write(lines,file="targetInhibitionNetwork.sif", sep = "\n", append=T)
    }
    
    # connecting the last component with T2
    for (j in expressions_compact[[num_components]]){
      lines <- paste(c(Terminals[2],"xx",j), collapse="\t")
      write(lines,file="targetInhibitionNetwork.sif", sep = "\n", append=T)
    }
    
    # connecting the middle components
    for (i in 1:(num_components-1)) {
      for (j in expressions_compact[[i]]){
        for (k in expressions_compact[[i+1]]){
          lines <- paste(c(j,"xx",k), collapse="\t")
          write(lines,file="targetInhibitionNetwork.sif", sep = "\n", append=T)
        }
      }
    } 
    cat()
    #dev.off()
} 
