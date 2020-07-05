# 3D plot using R rgl package
# adopted from http://www.sthda.com/english/wiki/a-complete-guide-to-3d-visualization-device-system-in-r-r-software-and-data-visualization
install.packages("rgl")

library("rgl")


rgl.open() # Open a new RGL device
rgl.bg(color = "white") # Setup the background color


# custom initialize function
#' @param new.device a logical value. If TRUE, creates a new device
#' @param bg the background color of the device
#' @param width the width of the device
rgl_init <- function(new.device = FALSE, bg = "white", width = 640) {
  if( new.device | rgl.cur() == 0 ) {
    rgl.open()
    par3d(windowRect = 50 + c( 0, 0, width, width ) )
    rgl.bg(color = bg )
  }
  rgl.clear(type = c("shapes", "bboxdeco"))
  rgl.viewpoint(theta = 15, phi = 20, zoom = 1)
}

#' Get colors for the different levels of
#' a factor variable
#'
#' @param groups a factor variable containing the groups
#'  of observations
#' @param colors a vector containing the names of
#   the default colors to be used
get_colors <- function(groups, group.col = palette()){
  groups <- as.factor(groups)
  ngrps <- length(levels(groups))
  if(ngrps > length(group.col))
    group.col <- rep(group.col, ngrps)
  color <- group.col[as.numeric(groups)]
  names(color) <- as.vector(groups)
  return(color)
}


# read in coordinates
xyz <- read.table('conformation1.txt')
# read in raw counts file to get row names
counts <- read.table('GG.count.matrix.noXY.csv', sep=',', head=T, row.names=1)
chrs <- rownames(counts)
rownames(xyz) <- chrs
# get the chr IDs
chrs <- paste('chr', substr(chrs[1:19], 2, 10), sep='')
# set the sphere size for each chr
sizes <- log2(colSums(counts))

# plot
rgl_init()
colors <- rainbow(19)
plot3d(xyz, xlab = '', ylab = '', zlab = '', type='s', radius=sizes, col=colors, aspect=T)
texts3d(xyz-20, text=names(sizes))
for (i in 1:19) {
	lines3d(xyz[c(i, i+19), ], col=colors[i])
}

rgl.postscript("plot.pdf",fmt="pdf")
