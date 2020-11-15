# Parameters:
# col: list of colours, e.g. c("#08306B","#2171B5","#6BAED6","#C6DBEF","#F7FBFF")
# brks: sequence of intervals, e.g. c(-0.90,-0.85,-0.80,-0.75,-0.70,-0.65)
# height: height of the rectangles, e.g. 1
# cex.axis: relative size of the label, default is 1

# Creating vertically aligned colorbar with custom colors:
colorbar_vertical <- function(col, brks, height, cex.axis=1) {
                      seq <- 1:length(brks)
                      plot(c(seq[1],tail(seq,1)), c(seq[1],tail(seq,1)),
                             pch="", type="n", bty="n", axes=FALSE, ann=FALSE) 
                        axis(2, at=seq, labels=brks, cex.axis=0.8, las=2)
                        for (i in 1:length(col)) {
                          rect(xleft=0, xright=height, ybottom=seq[i], ytop=seq[i+1], 
                               col=col[i], border=NA)
                        }
}

# Creating horizontally aligned colorscale with custom colors:
colorbar_horizontal <- function(col, brks, height, cex.axis=1) {
                      seq <- 1:length(brks)
                      plot(c(seq[1],tail(seq,1)), c(seq[1],tail(seq,1)),
                           pch="", type="n", bty="n", axes=FALSE, ann=FALSE) 
                        axis(1, at=seq, labels=brks, cex.axis=0.8)
                        for (i in 1:length(col)) {
                        rect(xleft=seq[i], xright=seq[i+1], ybottom=0, ytop=height, 
                             col=col[i], border=NA)
                        }
}
