variofit2 <- function (vario, ini.cov.pars, cov.model, fix.nugget = FALSE, 
                       nugget = 0, fix.kappa = TRUE, kappa = 0.5, simul.number = NULL, 
                       max.dist = vario$max.dist, weights, minimisation.function, 
                       limits = pars.limits(), messages, ...) 
{
    call.fc <- match.call()
    if (missing(messages)) 
        messages.screen <- as.logical(ifelse(is.null(getOption("geoR.messages")), 
                                             TRUE, getOption("geoR.messages")))
    else messages.screen <- messages
    if (length(class(vario)) == 0 || all(class(vario) != "variogram")) 
        warning("object vario should preferably be of the geoR's class \"variogram\"")
    if (!missing(ini.cov.pars)) {
        if (any(class(ini.cov.pars) == "eyefit")) 
            cov.model <- ini.cov.pars[[1]]$cov.model
        if (any(class(ini.cov.pars) == "variomodel")) 
            cov.model <- ini.cov.pars$cov.model
    }
    if (missing(cov.model)) 
        cov.model <- "matern"
    cov.model <- match.arg(cov.model, choices = .geoR.cov.models)
    if (cov.model == "stable") 
        cov.model <- "powered.exponential"
    if (cov.model == "powered.exponential") 
        if (limits$kappa["upper"] > 2) 
            limits$kappa["upper"] <- 2
    if (missing(weights)) {
        if (vario$output.type == "cloud") 
            weights <- "equal"
        else weights <- "npairs"
    }
    else weights <- match.arg(weights, choices = c("npairs", 
                                                   "equal", "cressie"))
    if (messages.screen) {
        cat(paste("variofit: covariance model used is", cov.model, 
                  "\n"))
        cat(paste("variofit: weights used:", weights, "\n"))
    }
    if (missing(minimisation.function)) 
        minimisation.function <- "optim"
    if (any(cov.model == c("linear", "power")) & minimisation.function == 
        "nls") {
        cat("warning: minimisation function nls can not be used with given cov.model.\n          changing for \"optim\".\n")
        minimisation.function <- "optim"
    }
    if (minimisation.function == "nls" & weights != "equal") {
        warning("variofit: minimisation function nls can only be used with weights=\"equal\".\n          changing for \"optim\".\n")
        minimisation.function <- "optim"
    }
    if (is.matrix(vario$v) & is.null(simul.number)) 
        stop("object in vario$v is a matrix. This function works for only 1 empirical variogram at once\n")
    if (!is.null(simul.number)) 
        vario$v <- vario$v[, simul.number]
    if (mode(max.dist) != "numeric" || length(max.dist) > 1) 
        stop("a single numerical value must be provided in the argument max.dist")
    if (max.dist == vario$max.dist) 
        XY <- list(u = vario$u, v = vario$v, n = vario$n)
    else XY <- list(u = vario$u[vario$u <= max.dist], v = vario$v[vario$u <= 
                                                                      max.dist], n = vario$n[vario$u <= max.dist])
    if (cov.model == "pure.nugget") {
        minimisation.function <- "not used"
        message <- "correlation function does not require numerical minimisation"
        if (weights == "equal") 
            lm.wei <- rep(1, length(XY$u))
        else lm.wei <- XY$n
        if (cov.model == "pure.nugget") {
            if (fix.nugget) {
                temp <- lm((XY$v - nugget) ~ 1, weights = lm.wei)
                cov.pars <- c(temp$coef, 0)
            }
            else {
                temp <- lm(XY$v ~ 1, weights = lm.wei)
                nugget <- temp$coef
                cov.pars <- c(0, 0)
            }
        }
        value <- sum((temp$residuals)^2)
    }
    else {
        if (messages.screen) 
            cat(paste("variofit: minimisation function used:", 
                      minimisation.function, "\n"))
        umax <- max(vario$u)
        vmax <- max(vario$v)
        if (missing(ini.cov.pars)) {
            ini.cov.pars <- as.matrix(expand.grid(c(vmax/2, 3 * 
                                                        vmax/4, vmax), seq(0, 0.8 * umax, len = 6)))
            if (!fix.nugget) 
                nugget <- unique(c(nugget, vmax/10, vmax/4, vmax/2))
            if (!fix.kappa) 
                kappa <- unique(c(kappa, 0.25, 0.5, 1, 1.5, 2))
            if (messages.screen) 
                warning("initial values not provided - running the default search")
        }
        else {
            if (any(class(ini.cov.pars) == "eyefit")) {
                init <- nugget <- kappa <- NULL
                for (i in 1:length(ini.cov.pars)) {
                    init <- drop(rbind(init, ini.cov.pars[[i]]$cov.pars))
                    nugget <- c(nugget, ini.cov.pars[[i]]$nugget)
                    if (cov.model == "gneiting.matern") 
                        kappa <- drop(rbind(kappa, ini.cov.pars[[i]]$kappa))
                    else kappa <- c(kappa, ini.cov.pars[[i]]$kappa)
                }
                ini.cov.pars <- init
            }
            if (any(class(ini.cov.pars) == "variomodel")) {
                nugget <- ini.cov.pars$nugget
                kappa <- ini.cov.pars$kappa
                ini.cov.pars <- ini.cov.pars$cov.pars
            }
        }
        if (is.matrix(ini.cov.pars) | is.data.frame(ini.cov.pars)) {
            ini.cov.pars <- as.matrix(ini.cov.pars)
            if (nrow(ini.cov.pars) == 1) 
                ini.cov.pars <- as.vector(ini.cov.pars)
            else {
                if (ncol(ini.cov.pars) != 2) 
                    stop("\nini.cov.pars must be a matrix or data.frame with 2 components: \ninitial values for sigmasq (partial sill) and phi (range parameter)\n")
            }
        }
        else if (length(ini.cov.pars) > 2) 
            stop("\nini.cov.pars must provide initial values for sigmasq and phi\n")
        if (is.matrix(ini.cov.pars) | (length(nugget) > 1) | 
            (length(kappa) > 1)) {
            if (messages.screen) 
                cat("variofit: searching for best initial value ...")
            ini.temp <- matrix(ini.cov.pars, ncol = 2)
            grid.ini <- as.matrix(expand.grid(sigmasq = unique(ini.temp[, 
                                                                        1]), phi = unique(ini.temp[, 2]), tausq = unique(nugget), 
                                              kappa = unique(kappa)))
            v.loss <- function(parms, u, v, n, cov.model, weights) {
                sigmasq <- parms[1]
                phi <- parms[2]
                if (cov.model == "power") 
                    phi <- 2 * exp(phi)/(1 + exp(phi))
                tausq <- parms[3]
                kappa <- parms[4]
                if (cov.model == "power") 
                    v.mod <- tausq + cov.spatial(u, cov.pars = c(sigmasq, 
                                                                 phi), cov.model = "power", kappa = kappa)
                else v.mod <- (sigmasq + tausq) - cov.spatial(u, 
                                                              cov.pars = c(sigmasq, phi), cov.model = cov.model, 
                                                              kappa = kappa)
                if (weights == "equal") 
                    loss <- sum((v - v.mod)^2)
                if (weights == "npairs") 
                    loss <- sum(n * (v - v.mod)^2)
                if (weights == "cressie") 
                    loss <- sum((n/(v.mod^2)) * (v - v.mod)^2)
                return(loss)
            }
            grid.loss <- apply(grid.ini, 1, v.loss, u = XY$u, 
                               v = XY$v, n = XY$n, cov.model = cov.model, weights = weights)
            ini.temp <- grid.ini[which(grid.loss == min(grid.loss))[1], 
                                 , drop = FALSE]
            if (is.R()) 
                rownames(ini.temp) <- "initial.value"
            if (messages.screen) {
                cat(" selected values:\n")
                print(rbind(round(ini.temp, digits = 2), status = ifelse(c(FALSE, 
                                                                           FALSE, fix.nugget, fix.kappa), "fix", "est")))
                cat(paste("loss value:", min(grid.loss), "\n"))
            }
            names(ini.temp) <- NULL
            ini.cov.pars <- ini.temp[1:2]
            nugget <- ini.temp[3]
            kappa <- ini.temp[4]
            grid.ini <- NULL
        }
        if (ini.cov.pars[1] > 2 * vmax) 
            warning("unreasonable initial value for sigmasq (too high)")
        if (ini.cov.pars[1] + nugget > 3 * vmax) 
            warning("unreasonable initial value for sigmasq + nugget (too high)")
        if (vario$output.type != "cloud") {
            if (ini.cov.pars[1] + nugget < 0.3 * vmax) 
                warning("unreasonable initial value for sigmasq + nugget (too low)")
        }
        if (nugget > 2 * vmax) 
            warning("unreasonable initial value for nugget (too high)")
        if (ini.cov.pars[2] > 1.5 * umax) 
            warning("unreasonable initial value for phi (too high)")
        if (!fix.kappa) {
            if (cov.model == "powered.exponential") 
                Tkappa.ini <- log(kappa/(2 - kappa))
            else Tkappa.ini <- log(kappa)
        }
        if (minimisation.function == "nls") {
            if (ini.cov.pars[2] == 0) 
                ini.cov.pars <- max(XY$u)/10
            if (kappa == 0) 
                kappa <- 0.5
            if (cov.model == "power") 
                Tphi.ini <- log(ini.cov.pars[2]/(2 - ini.cov.pars[2]))
            else Tphi.ini <- log(ini.cov.pars[2])
            XY$cov.model <- cov.model
            if (fix.nugget) {
                XY$nugget <- as.vector(nugget)
                if (fix.kappa) {
                    XY$kappa <- as.vector(kappa)
                    res <- nls((v - nugget) ~ matrix((1 - cov.spatial(u, 
                                                                      cov.pars = c(1, exp(Tphi)), cov.model = cov.model, 
                                                                      kappa = kappa)), ncol = 1), start = list(Tphi = Tphi.ini), 
                               data = XY, algorithm = "plinear", ...)
                }
                else {
                    if (cov.model == "powered.exponential") 
                        res <- nls((v - nugget) ~ matrix((1 - cov.spatial(u, 
                                                                          cov.pars = c(1, exp(Tphi)), cov.model = cov.model, 
                                                                          kappa = (2 * exp(Tkappa)/(1 + exp(Tkappa))))), 
                                                         ncol = 1), start = list(Tphi = Tphi.ini, 
                                                                                 Tkappa = Tkappa.ini), data = XY, algorithm = "plinear", 
                                   ...)
                    else res <- nls((v - nugget) ~ matrix((1 - 
                                                               cov.spatial(u, cov.pars = c(1, exp(Tphi)), 
                                                                           cov.model = cov.model, kappa = exp(Tkappa))), 
                                                          ncol = 1), start = list(Tphi = Tphi.ini, 
                                                                                  Tkappa = Tkappa.ini), data = XY, algorithm = "plinear", 
                                    ...)
                    kappa <- exp(coef(res)["Tkappa"])
                    names(kappa) <- NULL
                }
                cov.pars <- coef(res)[c(".lin", "Tphi")]
                names(cov.pars) <- NULL
            }
            else {
                if (fix.kappa) {
                    XY$kappa <- kappa
                    res <- nls(v ~ cbind(1, (1 - cov.spatial(u, 
                                                             cov.pars = c(1, exp(Tphi)), cov.model = cov.model, 
                                                             kappa = kappa))), start = list(Tphi = Tphi.ini), 
                               algorithm = "plinear", data = XY, ...)
                }
                else {
                    if (cov.model == "powered.exponential") 
                        res <- nls(v ~ cbind(1, (1 - cov.spatial(u, 
                                                                 cov.pars = c(1, exp(Tphi)), cov.model = cov.model, 
                                                                 kappa = (2 * exp(Tkappa)/(1 + exp(Tkappa)))))), 
                                   start = list(Tphi = Tphi.ini, Tkappa = Tkappa.ini), 
                                   algorithm = "plinear", data = XY, ...)
                    else res <- nls(v ~ cbind(1, (1 - cov.spatial(u, 
                                                                  cov.pars = c(1, exp(Tphi)), cov.model = cov.model, 
                                                                  kappa = exp(Tkappa)))), start = list(Tphi = Tphi.ini, 
                                                                                                       Tkappa = Tkappa.ini), algorithm = "plinear", 
                                    data = XY, ...)
                    kappa <- exp(coef(res)["Tkappa"])
                    names(kappa) <- NULL
                }
                nugget <- coef(res)[".lin1"]
                names(nugget) <- NULL
                cov.pars <- coef(res)[c(".lin2", "Tphi")]
                names(cov.pars) <- NULL
            }
            if (cov.model == "power") 
                cov.pars[2] <- 2 * exp(cov.pars[2])/(1 + exp(cov.pars[2]))
            else cov.pars[2] <- exp(cov.pars[2])
            if (nugget < 0 | cov.pars[1] < 0) {
                warning("\nvariofit: negative variance parameter found using the default option \"nls\".\n        Try another minimisation function and/or fix some of the parameters.\n")
                temp <- c(sigmasq = cov.pars[1], phi = cov.pars[2], 
                          tausq = nugget, kappa = kappa)
                print(rbind(round(temp, digits = 4), status = ifelse(c(FALSE, 
                                                                       FALSE, fix.nugget, fix.kappa), "fix", "est")))
                return(invisible())
            }
            value <- sum(resid(res)^2)
            message <- "nls does not provides convergence message"
        }
        if (minimisation.function == "nlm" | minimisation.function == 
            "optim") {
            .global.list <- list(u = XY$u, v = XY$v, n = XY$n, 
                                 fix.nugget = fix.nugget, nugget = nugget, fix.kappa = fix.kappa, 
                                 kappa = kappa, cov.model = cov.model, m.f = minimisation.function, 
                                 weights = weights)
            ini <- ini.cov.pars
            if (cov.model == "power") 
                ini[2] <- log(ini[2]/(2 - ini[2]))
            if (cov.model == "linear") 
                ini <- ini[1]
            if (fix.nugget == FALSE) 
                ini <- c(ini, nugget)
            if (!fix.kappa) 
                ini <- c(ini, Tkappa.ini)
            names(ini) <- NULL
            if (minimisation.function == "nlm") {
                result <- nlm(.loss.vario, ini, g.l = .global.list, 
                              ...)
                result$par <- result$estimate
                result$value <- result$minimum
                result$convergence <- result$code
                if (!is.null(get(".temp.theta", pos = 1))) 
                    result$par <- get(".temp.theta", pos = 1)
            }
            else {
                lower.l <- sapply(limits, function(x) x[1])
                upper.l <- sapply(limits, function(x) x[2])
                if (fix.kappa == FALSE) {
                    if (fix.nugget) {
                        lower <- lower.l[c("sigmasq.lower", "phi.lower", 
                                           "kappa.lower")]
                        upper <- upper.l[c("sigmasq.upper", "phi.upper", 
                                           "kappa.upper")]
                    }
                    else {
                        lower <- lower.l[c("sigmasq.lower", "phi.lower", 
                                           "tausq.rel.lower", "kappa.lower")]
                        upper <- upper.l[c("sigmasq.upper", "phi.upper", 
                                           "tausq.rel.upper", "kappa.upper")]
                    }
                }
                else {
                    if (cov.model == "power") {
                        if (fix.nugget) {
                            lower <- lower.l[c("sigmasq.lower", "phi.lower")]
                            upper <- upper.l[c("sigmasq.upper", "phi.upper")]
                        }
                        else {
                            lower <- lower.l[c("sigmasq.lower", "phi.lower", 
                                               "tausq.rel.lower")]
                            upper <- upper.l[c("sigmasq.upper", "phi.upper", 
                                               "tausq.rel.upper")]
                        }
                    }
                    else {
                        lower <- lower.l["phi.lower"]
                        upper <- upper.l["phi.upper"]
                    }
                }
                result <- optim(ini, .loss.vario, method = "L-BFGS-B", 
                                hessian = TRUE, lower = lower, upper = upper, 
                                g.l = .global.list, ...)
            }
            value <- result$value
            message <- paste(minimisation.function, "convergence code:", 
                             result$convergence)
            if (cov.model == "linear") 
                result$par <- c(result$par[1], 1, result$par[-1])
            cov.pars <- as.vector(result$par[1:2])
            if (cov.model == "power") 
                cov.pars[2] <- 2 * exp(cov.pars[2])/(1 + exp(cov.pars[2]))
            if (!fix.kappa) {
                if (fix.nugget) 
                    kappa <- result$par[3]
                else {
                    nugget <- result$par[3]
                    kappa <- result$par[4]
                }
                if (.global.list$cov.model == "powered.exponential") 
                    kappa <- 2 * (exp(kappa))/(1 + exp(kappa))
                else kappa <- exp(kappa)
            }
            else if (!fix.nugget) 
                nugget <- result$par[3]
        }
    }
    estimation <- list(nugget = nugget, cov.pars = cov.pars, 
                       cov.model = cov.model, kappa = kappa, value = value, 
                       trend = vario$trend, beta.ols = vario$beta.ols, practicalRange = practicalRange(cov.model = cov.model, 
                                                                                                       phi = cov.pars[2], kappa = kappa), max.dist = max.dist, 
                       minimisation.function = minimisation.function)
    estimation$weights <- weights
    if (weights == "equal") 
        estimation$method <- "OLS"
    else estimation$method <- "WLS"
    estimation$fix.nugget <- fix.nugget
    estimation$fix.kappa <- fix.kappa
    estimation$lambda <- vario$lambda
    estimation$message <- message
    estimation$call <- call.fc
    estimation$result <- result
    oldClass(estimation) <- c("variomodel", "variofit")
    estimation$hessian <- result$hessian
    return(estimation)
}
