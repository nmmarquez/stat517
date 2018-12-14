rm(list=ls())
library(INLA)
library(INSP)
library(dplyr)
library(rgeos)

mesh <- gCentroid(mx.sp.df, byid=T) %>% inla.mesh.create
alpha = 2 # refers to v = a - d/2
param = NULL 
constr = FALSE 
extraconstr.int = NULL 
extraconstr = NULL
fractional.method = c("parsimonious", "null") # lets just deal with integers for now
B.tau = matrix(c(0, 1, 0), 1, 3) # the log linear model for tau 
B.kappa = matrix(c(0, 0, 1), 1, 3) # log linear model for kappa
prior.variance.nominal = 1 
prior.range.nominal = NULL 
prior.tau = NULL
prior.kappa = NULL 
theta.prior.mean = NULL 
theta.prior.prec = 0.1 

fractional.method = "null"
if (!is.null(param)){
    deprecated = !c(missing(B.tau), missing(B.kappa), missing(prior.variance.nominal), 
                    missing(prior.range.nominal), missing(prior.tau), 
                    missing(prior.kappa), missing(theta.prior.mean), 
                    missing(theta.prior.prec))
    deprecated = c("B.tau", "B.kappa", "prior.variance.nominal", 
                   "prior.range.nominal", "prior.tau", "prior.kappa", 
                   "theta.prior.mean", "theta.prior.prec")[deprecated]
    if (length(deprecated) > 0) {
        warning(paste("'param' specified;  ", "Ignoring deprecated parameter(s) ", 
                      paste(deprecated, collapse = ", "), ".", sep = ""))
    }
}
if (is.null(param)) {
    param = param2.matern.orig(mesh, alpha, B.tau, B.kappa, 
                               prior.variance.nominal, prior.range.nominal, prior.tau, 
                               prior.kappa, theta.prior.mean, theta.prior.prec)
}

# detrmine the dimensionality of the mesh
d = ifelse(inherits(mesh, "inla.mesh"), 2, 1)
# calculate nu as specified in equation 2
nu = alpha - d/2
# some back up values in case nu is less than .5
nu.nominal = max(0.5, nu)
alpha.nominal = max(nu.nominal + d/2, alpha)
n.spde = ifelse(d == 2, mesh$n, mesh$m)
n.theta = ncol(B.kappa) - 1L
# go through this sequence if the dimensionality is 2
if (d == 2) {
    fem = inla.fmesher.smorg(mesh$loc, mesh$graph$tv, fem = 2,
                             output = list("c0", "c1", "g1", "g2"))
    # loc = mesh$loc
    # tv = mesh$graph$tv
    # fem = NULL
    # aniso = NULL
    # gradients = FALSE
    # sph0 = NULL
    # sph = NULL
    # bspline = NULL
    # points2mesh = NULL
    # splitlines = NULL
    # output = list("c0", "c1", "g1", "g2")
    # keep = FALSE
    # prefix = INLA:::inla.fmesher.make.prefix(NULL, NULL)
    # n = nrow(loc)
    # s.dim = ncol(loc)
    # output.given = TRUE
    # output.fem = list("c0", "g1", "g2")
    # output.aniso = list("g1aniso", "g2aniso")
    # output.gradients = list("dx", "dy", "dz")
    # output.sph0 = list("sph0")
    # output.sph = list("sph")
    # output.bspline = list("bspline")
    # output.p2m = list("p2m.t", "p2m.b")
    # output.splitlines <- list("split.loc", "split.idx", "split.origin", 
    #                           "split.t", "split.b1", "split.b2")
    # indexoutput <- list("split.idx", "split.t", "split.origin")
    # INLA:::fmesher.write(INLA:::inla.affirm.double(loc), prefix, "s")
    # INLA:::fmesher.write(INLA:::inla.affirm.integer(tv) - 1L, prefix, "tv")
    # all.args = "--smorg --input=s,tv"
    # all.args = paste(all.args, inla.getOption("fmesher.arg"))
    # echoc = INLA:::inla.fmesher.call(all.args = all.args, prefix = prefix)
    # result = list()
    # for (name in output) {
    #     if (identical(name, "p2m.t")) 
    #         if (!file.exists(paste(prefix, name, sep = ""))) 
    #             result[[name]] = INLA:::fmesher.read(prefix, "points2mesh.t") + 
    #                 1L
    #         else result[[name]] = INLA:::fmesher.read(prefix, name) + 
    #                 1L
    #     else if (identical(name, "p2m.b")){
    #         if (!file.exists(paste(prefix, name, sep = ""))){
    #             result[[name]] = INLA:::fmesher.read(prefix, "points2mesh.b")
    #         }
    #         else{ 
    #             result[[name]] = INLA:::fmesher.read(prefix, name)
    #         }
    #     }
    #     else if(name %in% indexoutput){
    #         result[[name]] = INLA:::fmesher.read(prefix, name) + 1L
    #     }
    #     else{
    #         result[[name]] = INLA:::fmesher.read(prefix, name)
    #     }
    # }
    # if (!keep){
    #     unlink(paste(prefix, "*", sep = ""), recursive = FALSE)
    # }
    # fem = result
}
# go through this sequence if the dimensionality is 1
if(d ==1) {
    fem = inla.mesh.1d.fem(mesh)
    if (mesh$degree == 2) {
        fem$c0 = fem$c1
    }
}
if (alpha == 2) {
    B.phi0 = param$B.tau
    B.phi1 = 2 * param$B.kappa
    M0 = fem$c0
    M1 = fem$g1
    M2 = fem$g2
}
if (alpha == 1) {
    B.phi0 = param$B.tau
    B.phi1 = param$B.kappa
    M0 = fem$c0
    M1 = fem$g1 * 0
    M2 = fem$g1
}


spde = inla.spde2.generic(M0 = M0, M1 = M1, M2 = M2, 
                          B0 = B.phi0, B1 = B.phi1, B2 = 1, theta.mu = param$theta.prior.mean, 
                          theta.Q = param$theta.prior.prec, transform = "identity", 
                          BLC = param$BLC)

B0_ = B.phi0
B1_ = B.phi1
B2_ = 1
theta.mu = param$theta.prior.mean
theta.Q = param$theta.prior.prec
transform = "identity"
BLC_ = param$BLC
theta.initial = theta.mu
fixed = rep(FALSE, length(theta.mu))
theta.fixed = theta.initial[fixed]

M0 = INLA:::inla.as.dgTMatrix(M0)
M1 = INLA:::inla.as.dgTMatrix(M1)
M2 = INLA:::inla.as.dgTMatrix(M2)
n.spde = nrow(M0)
n.theta = length(theta.mu)
spde = list(
    model = "generic", 
    n.spde = n.spde, 
    n.theta = n.theta, 
    param.inla = list(), 
    f = list())
class(spde) = c("inla.spde2", "inla.spde", "inla.model.class")
B0 = INLA:::inla.spde.homogenise_B_matrix(B0_, n.spde, n.theta)
B1 = INLA:::inla.spde.homogenise_B_matrix(B1_, n.spde, n.theta)
B2 = INLA:::inla.spde.homogenise_B_matrix(B2_, n.spde, n.theta)
BLC = INLA:::inla.spde.homogenise_B_matrix(BLC_, nrow(BLC_), n.theta)
param.inla = list(n = n.spde, n.theta = n.theta, M0 = M0, 
                  M1 = M1, M2 = M2, B0 = B0, B1 = B1, B2 = B2, BLC = BLC, 
                  theta.mu = theta.mu, theta.Q = theta.Q, transform = transform, 
                  theta.initial = theta.initial, fixed = fixed, theta.fixed = theta.fixed)

param.inla$BLC = param.inla$BLC[rowSums(abs(param.inla$BLC[, -1, drop = FALSE])) > 0, , drop = FALSE]

spde$param.inla = param.inla
spde$f = (list(model = "spde2", n = n.spde, spde2.transform = transform, 
               hyper.default = (list(theta1 = list(prior = "mvnorm", 
                                                   param = (c(param.inla$theta.mu, as.matrix(param.inla$theta.Q))))))))
 
for (k in 1:n.theta) {
    eval(parse(text = paste("spde$f$hyper.default$theta", 
                            k, "$initial = param.inla$theta.mu[k]", sep = "")))
}


spde$model = "matern"
spde$BLC = param$BLC

spde$mesh <- mesh
invisible(spde)

