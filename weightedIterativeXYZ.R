countsFile <- 'G1G2.count.noY.normed.txt'

# read in interaction frequency file
counts <- read.table(countsFile, head=T, row.names=1, sep='\t')
ncols <- ncol(counts)
counts <- counts[, -ncols] # remove genotype
nrows <- nrow(counts)
counts <- as.matrix(counts)
diag(counts) <- 0


# mask out homologous chromosomes
half <- nrows / 2
for (i in 1:nrows) {
    if (i > half) {
        counts[i, i - half] <- 0
    }
    else {
        counts[i, i + half] <- 0
    }
}

# assume homologous chromosomes have the longest/mean distance
half <- nrows / 2
for (i in 1:nrows) {
    if (i > half) {
        # counts[i, i - half] <- sort(counts[i, ])[2]
        counts[i, i - half] <- mean(counts[i, counts[i, ] > 0])
    }
    else {
        # counts[i, i + half] <- sort(counts[i, ])[2]
        counts[i, i + half] <- mean(counts[i, counts[i, ] > 0])
    }
}

normed.log2 <- log2(counts + 1)

cSum <- colSums(counts)
allchrs <- rownames(counts)

g1 <- read.table('G1.dis2count.txt')
g1 <- g1[-1,]
g1[,3] <- log2(g1[,2])
fg1 <- lowess(g1[,c(1,3)], f=0.01)

# calculate distance
distance <- t(apply(normed.log2, 1, function(x) {ind <- c(); for (y in x){ind <- c(ind, which(fg1$y < y)[1])}; return(ind)}))
rownames(distance) <- allchrs
colnames(distance) <- allchrs
diag(distance) <- 0


# mask out homologous chromosomes
# for (i in 1:nrows) {
#     if (i > half) {
#         distance[i, i - half] <- 0
#     }
#     else {
#         distance[i, i + half] <- 0
#     }
# }


distanceSqs <- distance ** 2
weights <- counts / rowSums(counts)


# initialize coordinates
xyz <- matrix(runif(3 * nrows, min=10, max=100), ncol=3, nrow=nrows)
xyz.init <- xyz

xyzDis <- apply(xyz, 1, function(x){rowSums(t(x - t(xyz)) ** 2)})
errDisSum <- sum(abs(distanceSqs - xyzDis) * weights)

# iterative correction
maxIter <- 2000
for (iter in 1:maxIter) {
	for (i in 1:nrows) {
		ratio <- distanceSqs / xyzDis
		diag(ratio) <- 0
		xyzFactor <- t(xyz[i, ] - t(xyz)) * ratio[i, ] + xyz
		xyz[i, ] <- colSums((t(xyz[i, ] - t(xyz)) * ratio[i, ] + xyz) * weights[i, ])
		xyzDis <- apply(xyz, 1, function(x){rowSums(t(x - t(xyz)) ** 2)})
		errDisSum <- c(sum(abs(distanceSqs - xyzDis) * weights), errDisSum)
		if (min(errDisSum) == errDisSum[1]) {
			xyz.best <- xyz
		}
	}
}
min(errDisSum) / sum(distanceSqs) * 100

xyz <- xyz.best
rownames(xyz) <- allchrs
colnames(xyz) <- c('x', 'y', 'z')
xyzDis <- apply(xyz, 1, function(x){rowSums(t(x - t(xyz)) ** 2)})


errDisSum <- rev(errDisSum) / sum(distanceSqs) * 100 # variance%
ld <- length(errDisSum)

pdf(paste('weighted.variance', nrows, maxIter, 'pdf', sep='.'))

x <- 1
plot(x=x:ld, y=errDisSum[x:ld], pch=20, xlab='Iteration', ylab='Variance(%)', main=paste(x, ':', ld, sep=''))
x <- 500
plot(x=x:ld, y=errDisSum[x:ld], pch=20, xlab='Iteration', ylab='Variance(%)', main=paste(x, ':', ld, sep=''))

dev.off()



# plot3d
library("rgl")
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
rgl_init()


half <- nrows / 2
colors <- rainbow(half)

# set the sphere size for each chr
sizes <- log10(cSum) * 10

plot3d(xyz, xlab='', ylab='', zlab='', type='s', radius=sizes, col=colors, aspect=T)

nlines <- 5
for (i in 1:nrows) {
	xyz.distance <- 1 / rowSums(t((xyz[i,] - t(xyz)) ** 2))
	ordered <- order(xyz.distance, decreasing=T)[1:nlines + 1]
	xyz.sum <- sum(xyz.distance[ordered])
	for (j in ordered) {
		d <- xyz.distance[j] / xyz.sum * 5 # enlarge the line width
		material3d(lwd=d)
		if (i > half) {
			lines3d(xyz[c(i, j), ], col=colors[i - half])
		}
		else {
			lines3d(xyz[c(i, j), ], col=colors[i])
		}
	}
}



rgl.postscript(paste("xyz", nrows, maxIter, "meanDis.3.pdf", sep='.'),fmt="pdf")

plot3d(xyz, xlab='', ylab='', zlab='', type='s', radius=sizes, col=colors, aspect=T)
texts3d(xyz+60, text=allchrs, cex=1)

rgl.postscript(paste("xyz", nrows, maxIter, "meanDis.4.pdf", sep='.'),fmt="pdf")

rownames(xyz.init) <- allchrs
colnames(xyz.init) <- c('x', 'y', 'z')
write.table(xyz.init, file=paste("xyz", nrows, maxIter, "meanDis.orig", sep='.'), quote=F, sep='\t')
write.table(xyz, file=paste("xyz", nrows, maxIter, 'meanDis.final', sep='.'), quote=F, sep='\t')
