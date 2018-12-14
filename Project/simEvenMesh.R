rm(list=ls())
library(dplyr)
library(Matrix)
library(ar.matrix) # devtools::install_github("nmmarquez/ar.matrix")
library(ggplot2)
library(INLA)

buildQv1 <- function(gridDim, kappa=.3, verbose=F, sparse=F){
    gridN <- gridDim^2
    matrixPos <- expand.grid(x=1:gridDim, y=1:gridDim) %>%
        mutate(Qpos=1:gridN)
    gridQ <- matrix(0, nrow=gridN, ncol=gridN)
    a_ <- kappa^2 + 4
    
    for(i in 1:gridN){
        if(verbose){
            print(i)
        }
        diagPos <- subset(matrixPos, Qpos==i)
        pos_ <- c(diagPos$x, diagPos$y)
        gridQ[i, i] <- 4 + a_^2
        for(j in list(c(1,0), c(0,1), c(-1,0), c(0,-1))){
            adj_pos <- pos_ + j
            if(all(adj_pos > 0) & all(adj_pos <= gridN)){
                Qpos2 <- subset(matrixPos, x==adj_pos[1] & y==adj_pos[2])$Qpos
                gridQ[i, Qpos2] <- -2 * a_
                gridQ[Qpos2, i] <- -2 * a_
            }
        }
        for(j in list(c(1,1), c(-1,1), c(-1,-1), c(1,-1))){
            adj_pos <- pos_ + j
            if(all(adj_pos > 0) & all(adj_pos <= gridN)){
                Qpos2 <- subset(matrixPos, x==adj_pos[1] & y==adj_pos[2])$Qpos
                gridQ[i, Qpos2] <- 2
                gridQ[Qpos2, i] <- 2 
            }
        }
        for(j in list(c(2,0), c(0,2), c(-2,0), c(0,-2))){
            adj_pos <- pos_ + j
            if(all(adj_pos > 0) & all(adj_pos <= gridN)){
                Qpos2 <- subset(matrixPos, x==adj_pos[1] & y==adj_pos[2])$Qpos
                gridQ[i, Qpos2] <- 1
                gridQ[Qpos2, i] <- 1
            }
        }
    }
    if(sparse){
        gridQ <- Matrix(gridQ, sparse=TRUE)
    }
    return(gridQ)
}

buildQv2 <- function(gridDim, kappa=.3, verbose=FALSE, sparse=FALSE){
    gridN <- gridDim^2
    matrixPos <- expand.grid(x=1:gridDim, y=1:gridDim) %>%
        mutate(Qpos=1:gridN)
    gridQ <- matrix(0, nrow=gridN, ncol=gridN)
    a_ <- kappa^2 + 4
    
    for(i in 1:gridN){
        if(verbose){
            print(i)
        }
        diagPos <- subset(matrixPos, Qpos==i)
        pos_ <- c(diagPos$x, diagPos$y)
        gridQ[i, i] <- a_ * (12 + a_^2)
        for(j in list(c(1,0), c(0,1), c(-1,0), c(0,-1))){
            adj_pos <- pos_ + j
            if(all(adj_pos > 0) & all(adj_pos <= gridN)){
                Qpos2 <- subset(matrixPos, x==adj_pos[1] & y==adj_pos[2])$Qpos
                gridQ[i, Qpos2] <- -3 * (3 + a_^2)
                gridQ[Qpos2, i] <- -3 * (3 + a_^2)
            }
        }
        for(j in list(c(1,1), c(-1,1), c(-1,-1), c(1,-1))){
            adj_pos <- pos_ + j
            if(all(adj_pos > 0) & all(adj_pos <= gridN)){
                Qpos2 <- subset(matrixPos, x==adj_pos[1] & y==adj_pos[2])$Qpos
                gridQ[i, Qpos2] <- 6 * a_
                gridQ[Qpos2, i] <- 6 * a_
            }
        }
        for(j in list(c(2,0), c(0,2), c(-2,0), c(0,-2))){
            adj_pos <- pos_ + j
            if(all(adj_pos > 0) & all(adj_pos <= gridN)){
                Qpos2 <- subset(matrixPos, x==adj_pos[1] & y==adj_pos[2])$Qpos
                gridQ[i, Qpos2] <- 3 * a_
                gridQ[Qpos2, i] <- 3 * a_
            }
        }
        for(j in list(c(3,0), c(0,3), c(-3,0), c(0,-3))){
            adj_pos <- pos_ + j
            if(all(adj_pos > 0) & all(adj_pos <= gridN)){
                Qpos2 <- subset(matrixPos, x==adj_pos[1] & y==adj_pos[2])$Qpos
                gridQ[i, Qpos2] <- -1
                gridQ[Qpos2, i] <- -1
            }
        }
        for(j in list(c(2,1), c(1,2), c(2,-1), c(-1,2), 
                      c(-2,1), c(1,-2), c(-2,-1), c(-1,-2))){
            adj_pos <- pos_ + j
            if(all(adj_pos > 0) & all(adj_pos <= gridN)){
                Qpos2 <- subset(matrixPos, x==adj_pos[1] & y==adj_pos[2])$Qpos
                gridQ[i, Qpos2] <- -3
                gridQ[Qpos2, i] <- -3
            }
        }
    }
    if(sparse){
        gridQ <- Matrix(gridQ, sparse=TRUE)
    }
    return(gridQ)
}

simQ <- function(Q){
    cholL <- chol(Q)
    z <- sapply(1, function(x) rnorm(nrow(Q)))
    x <- solve(t(cholL), z)@x
    return(x)
}

maternCor <- function(distX, kappa_, nu_=1){
    t3 <- besselK(kappa_ * distX, nu_)
    t2 <- (kappa_ * distX)^nu_
    t1 <- (2)^(1-nu_)/gamma(nu_)
    return(t1 * t2 * t3)
}

createDistGrid <- function(gridDim){
    gridN <- gridDim^2
    matrixPos <- expand.grid(x=1:gridDim, y=1:gridDim) %>%
        mutate(Qpos=1:gridN)
    gridSigma <- matrix(0, nrow=gridN, ncol=gridN)
    for(i in 1:gridN){
        for(j in i:gridN){
            p1 <- subset(matrixPos, Qpos == i)
            p2 <- subset(matrixPos, Qpos == j)
            dist_ <- c(dist(rbind(c(p1$x, p1$y), c(p2$x, p2$y))))
            gridSigma[i,j] <- dist_
            gridSigma[j,i] <- dist_
        }
    }
    return(gridSigma)
}

buildSigma <- function(gridDim, kappa=.3, nu_=1){
    gridSigma <- createDistGrid(gridDim)
    gridSigma <- maternCor(gridSigma, kappa, nu_=nu_)
    diag(gridSigma) <- rep(1, gridDim^2)
    return(gridSigma)
}

image(buildSigma(10, kappa=.9), 
      main="Exact Matern Covariance(nu=1)")
image(solve(buildQv1(10, kappa=.9)), 
      main="GMRF approximation of Matern Covariance(nu=1)")
image(solve(buildSigma(10, kappa=.9)), 
      main="Exact Matern Precision(nu=1)")
image(buildQv1(10, kappa=.9), 
      main="GMRF approximation of Matern Precision(nu=1)")
image(solve(buildSigma(10, nu_=2, kappa=.9)), 
      main="Exact Matern Precision(nu=2)")
image(buildQv2(10, kappa=.9), 
      main="GMRF approximation of Matern Precision(nu=2)")

rpos <- expand.grid(x=1:15, y=1:15) %>% mutate(pos=1:(15^2)) %>%
    filter(x%in%c(1,2,14,15) | y%in%c(1,2,14,15)) %>% select(pos) %>%
    unlist

m <- 60
system.time(Q1 <- buildQv1(m, kappa=.1, sparse=T))
system.time(Q2 <- buildQv2(m, kappa=.1, sparse=T))

set.seed(123)
expand.grid(x=1:m, y=1:m) %>%
    mutate(val_=rnorm(m^2)) %>%
    ggplot(aes(x=x, y=y, color=val_)) +
    scale_color_distiller(palette = "Spectral") +
    geom_point(size=2.5) + 
    labs(x="", y="", title="Points Simulated at Random") +
    theme_void()

expand.grid(x=1:m, y=1:m) %>%
    mutate(val_=simQ(Q1)) %>%
    ggplot(aes(x=x, y=y, color=val_)) +
    scale_color_distiller(palette = "Spectral") +
    geom_point(size=2.5) + 
    labs(x="", y="", title="Points Simulated with v=1") +
    theme_void()

expand.grid(x=1:m, y=1:m) %>%
    mutate(val_=simQ(Q2)) %>%
    ggplot(aes(x=x, y=y, color=val_)) +
    scale_color_distiller(palette = "Spectral") +
    geom_point(size=2.5) + 
    labs(x="", y="", title="Points Simulated with v=2") +
    theme_void()


pointsDF <- data.frame(x=runif(m^2,0, m), y=runif(m^2,0, m))

buildIrregQ <- function(pointsDF, kappa_=.3, nu_=2, sparse=TRUE){
    mesh <-  pointsDF %>% 
        as.matrix %>%
        inla.mesh.create
    fmesh <- inla.fmesher.smorg(mesh$loc, mesh$graph$tv, fem = 2,
                                output = list("c0", "c1", "g1", "g2"))
    M0 <- fmesh$c0
    M1 <- fmesh$g1
    M2 <- fmesh$g2
    Q <- as.matrix((kappa_^4)*M0 + 2*(kappa_^2) * M1 + M2)
    if(sparse){
        Q <- Matrix(Q, sparse=TRUE)
    }
    return(Q)
}

irregQ <- buildIrregQ(pointsDF, kappa_ = .1, sparse=T)
simX <- c(sim.AR(1, irregQ))

pointsDF %>%
    mutate(val_=simX[1:(m^2)]) %>%
    ggplot(aes(x=x, y=y, color=val_)) +
    scale_color_distiller(palette = "Spectral") +
    geom_point(size=2) + 
    labs(x="", y="", title="Points Simulated with v=1") +
    theme_void()

# use linear interpolation to project onto the full space
proj <- pointsDF %>%
    as.matrix %>%
    inla.mesh.create %>%
    inla.mesh.projector(dims=c(500,500))

pointsDF %>%
    as.matrix %>%
    inla.mesh.create %>%
    plot

data.frame(x=round(proj$lattice$loc[,1], 4), 
           y=round(proj$lattice$loc[,2], 4)) %>%
    mutate(z=c(inla.mesh.project(proj, field=simX))) %>%
    ggplot(aes(x, y, z=z)) + geom_tile(aes(fill = z)) + theme_bw() + 
    lims(y=c(0,m), x=c(0,m)) + 
    scale_fill_distiller(palette = "Spectral") +
    theme_void() + labs(title="Triangulated Interpolation")
