# colLab function of Joris Meys (http://stackoverflow.com/) modified by Ambrosio Torres to coloring edgePar and labels of a dendrogram 

colLab <- function(n)
{
  if(is.leaf(n)) {
    a <- attributes(n)
    # clusMember - a vector designating leaf grouping
    # labelColors - a vector of colors for the above grouping
    labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
    attr(n, "edgePar") <- list(col = labCol)
    attr(n, "label") <- NULL ;
  }
  n
}
