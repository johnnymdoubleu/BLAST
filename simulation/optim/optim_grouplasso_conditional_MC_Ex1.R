library(JOPS)
library(Pareto)
library(npreg)
library(gridExtra)
library(colorspace)
library(reshape2)
library(splines2)
library(scales)
suppressMessages(library(tidyverse))
# library(ggplotify)

#Scenario 1
# set.seed(12338)
# n <- 1000
# # beta <- c(0.2, 0.7)
# beta <- c(0.2, 0, 0.8, 0, 0, -0.1, 0, 0, 0, -0.4)
# p <- length(beta)
# x.scale <- cbind(replicate(p, runif(n, 0, 1)))
# x.scale <- scale(x.scale)


# alp.origin <- y <- NULL
# for(i in 1:n){
#   alp.origin[i] <- exp(sum(x.scale[i, ] * beta))
#   y[i] <- rPareto(1, 1, alpha = alp.origin[i])
# }
# # plot(y)

threshold <- 0.90
# x.origin <- x.scale
# y.origin <- y
# x.scale <- x.scale[which(y>quantile(y, threshold)),]
# u <- quantile(y, threshold)

# y <- as.matrix(y[y > quantile(y, threshold)])
# # n <- length(y)
# n <- length(y.origin)

# xholder <- new.bs.x <- bs.x <- matrix(,nrow=n, ncol=0)
# psi <- 33
# # new.bs.x <- bs.x <- array(NA, c(n, psi, p-1))
# for(i in 1:p){
#   splines <- bbase(x.origin[, i], min(x.origin[, i]), max(x.origin[, i]), nseg = (psi-3), bdeg = 3)
# #   splines <- bbase(x.scale[, i], min(x.scale[, i]), max(x.scale[, i]), nseg = (psi-3), bdeg = 3)
#   bs.x <- cbind(bs.x, splines)
# #   bs.x[, , i-1] <- splines
#   # phi <- dim(out[[1]][[1]]$X)[2]
#   # psi <- dim(splines)[2]
# }

# newx <- seq(0, 1, length.out=n)
# xholder <- cbind(xholder, rep(1,n))
# for(i in 1:p){
#   splines <- bbase(newx, min(newx), max(newx), nseg = (psi  -3), bdeg = 3)
#   new.bs.x <- cbind(new.bs.x, splines)
#   xholder <- cbind(xholder, newx)
# }
n <- 5000
psi <- 20

p <- 10
no.theta <- 2
simul.no <- 50


theta.container <- as.data.frame(matrix(, nrow = (no.theta *p), ncol= simul.no))
gamma.container <- as.data.frame(matrix(, nrow = (psi * p), ncol = simul.no))
# linear.container <- nonlinear.container <- f.container <- lapply(1:simul.no, matrix, data= NA, nrow=(n*(1-threshold)), ncol=p)
linear.container <- nonlinear.container <- f.container <- lapply(1:simul.no, data.frame)
alpha.container <- as.data.frame(matrix(, nrow=n, ncol = simul.no))

for (mm in 1:simul.no){
  #Scenario 1
    n <- 5000
    
    xholder.nonlinear <- xholder.linear <- bs.nonlinear <- bs.linear <- matrix(,nrow=n, ncol=0)
    x.origin <- cbind(replicate(p, runif(n, 0, 1)))
    for(i in 1:p){
        knots <- seq(min(x.origin[,i]), max(x.origin[,i]), length.out = psi)  
        tps <- basis.tps(x.origin[,i], knots, m = 2, rk = FALSE, intercept = TRUE)
        # tps <- mSpline(x.origin[,i], df=psi, Boundary.knots = range(x.origin[,i]), degree = 3, intercept=TRUE)
        #   bs.x <- cbind(bs.x, tps)
        bs.linear <- cbind(bs.linear, tps[,1:no.theta])
        bs.nonlinear <- cbind(bs.nonlinear, tps[,-c(1:no.theta)])  
    }

    gamma.origin <- matrix(, nrow = psi, ncol = p)
    for(j in 1:p){
        for (ps in 1:psi){
            if(j %in% c(2,4,5,6,9,10)){gamma.origin[ps, j] <- 0}
            else if(j==7){
                if(ps <= (psi/2)){gamma.origin[ps, j] <- 1}
                else{gamma.origin[ps, j] <- 1}
            }
            else {
                if(ps <= (psi/2)){gamma.origin[ps, j] <- 1}
                else{gamma.origin[ps, j] <- 1}
            }
        }
    }

    theta.origin <- matrix(, nrow = no.theta, ncol = p)
    for(j in 1:p){
        for (k in 1:no.theta){
            if(j %in% c(2,4,5,6,9,10)){theta.origin[k, j] <- 0}
            else if(j==7){
                if(k==1){theta.origin[k, j] <- 0.5}
                else{theta.origin[k, j] <- -0.3}
            }
            else {
                if(k==1){theta.origin[k,j] <- -0.2}
                else{theta.origin[k,j] <- 0.8}
            }
        }
    }

    f.nonlinear.origin <- f.linear.origin <- f.origin <- matrix(, nrow = n, ncol = p)
    for(j in 1:p){
        f.origin[, j] <- as.matrix(bs.linear[1:n, (((j-1)*no.theta)+1):(((j-1)*no.theta)+no.theta)]) %*% theta.origin[,j] + (bs.nonlinear[1:n,(((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j])
        f.linear.origin[,j] <- as.matrix(bs.linear[1:n, (((j-1)*no.theta)+1):(((j-1)*no.theta)+no.theta)]) %*% theta.origin[,j]
        f.nonlinear.origin[,j] <- (bs.nonlinear[1:n,(((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j])
    }

    alp.origin <- y.origin <- NULL
    for(i in 1:n){
        alp.origin[i] <- exp(sum(f.origin[i,]))
        y.origin[i] <- rPareto(1, 1, alpha = alp.origin[i])
    }

    u <- quantile(y.origin, threshold)
    x.origin <- x.origin[which(y.origin>u),]
    y.origin <- y.origin[y.origin > u]
    n <- length(y.origin)

    xholder.nonlinear <- xholder.linear <- bs.nonlinear <- bs.linear <- matrix(,nrow=n, ncol=0)
    newx <- seq(0, 1, length.out=n)
    xholder <- bs.x <- matrix(, nrow = n, ncol = p)
    for(i in 1:p){
        xholder[,i] <- seq(0, 1, length.out = n)
        test.knot <- seq(0, 1, length.out = psi)
        splines <- basis.tps(newx, test.knot, m=2, rk=FALSE, intercept = TRUE)
        xholder.linear <- cbind(xholder.linear, splines[,1:no.theta])
        xholder.nonlinear <- cbind(xholder.nonlinear, splines[,-c(1:no.theta)])
        knots <- seq(min(x.origin[,i]), max(x.origin[,i]), length.out = psi)  
        tps <- basis.tps(x.origin[,i], knots, m = 2, rk = FALSE, intercept = TRUE)
        # tps <- mSpline(x.origin[,i], df=psi, Boundary.knots = range(x.origin[,i]), degree = 3, intercept=TRUE)
        #   bs.x <- cbind(bs.x, tps)
        bs.linear <- cbind(bs.linear, tps[,1:no.theta])
        bs.nonlinear <- cbind(bs.nonlinear, tps[,-c(1:no.theta)])  
    }

    f.nonlinear.origin <- f.linear.origin <- f.origin <- matrix(, nrow = n, ncol = p)
    for(j in 1:p){
        f.origin[, j] <- as.matrix(bs.linear[1:n, (((j-1)*no.theta)+1):(((j-1)*no.theta)+no.theta)]) %*% theta.origin[,j] + (bs.nonlinear[1:n,(((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j])
        f.linear.origin[,j] <- as.matrix(bs.linear[1:n, (((j-1)*no.theta)+1):(((j-1)*no.theta)+no.theta)]) %*% theta.origin[,j]
        f.nonlinear.origin[,j] <- (bs.nonlinear[1:n,(((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j])
    }

    alp.origin <- NULL
    for(i in 1:n){
    alp.origin[i] <- exp(sum(f.origin[i,]))
    }

    # lambda <- 0.995#sqrt(2*log(n))
    # theta <- 0
    # lambda.1 <- lambda * theta
    # lambda.2 <-. lambda * (1-theta)
    lambda.1 <- 0.001
    lambda.2 <- 0
    lambda.3 <- 0.001

    log.posterior <- function(beta, y.origin){
        log.lik <- function(beta){
            exp.prime <- function(x, thres){
                if(x > thres){ans <- exp(thres) + exp(thres)*(x-thres)}
                else{ans <- exp(x)}
                return(ans)
            }
            theta <- matrix(beta[1:(no.theta*p)], nrow=no.theta)
            gamma <- matrix(beta[(no.theta*p)+1:(psi*p)], ncol=p)
            f <- matrix(, nrow=n, ncol=p)
            term <- first.term <- second.term <- NULL
            for(j in 1:p){
                # coef <- as.matrix(gamma[(((j-1)*psi)+1):(((j-1)*psi)+psi)])
                linear.term <- as.matrix(bs.linear[,(((j-1)*no.theta)+1):(((j-1)*no.theta)+no.theta)]) %*% theta[,j]
                nonlinear.term <- bs.nonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma[,j]
                f[, j] <- linear.term + nonlinear.term
            }
            for(i in 1:n){
                first.term[i] <- sum(f[i,]) - log(y.origin[i])
                second.term[i] <- exp.prime(sum(f[i,]), thres = 10) * log(y.origin[i]/u)
                term[i] <- first.term[i] - second.term[i]
            }
            return(sum(term))
        }
        log.prior <- function(beta){
            moreau.envelope <- function(w){
                if(w < -1){ans <- -0.5 - w}
                else if (1 < w){ans <- w - 0.5}
                else {ans <- (w^2) / 2}
                return(ans)
            }
            theta <- matrix(beta[1:(no.theta*p)], nrow=no.theta)
            gamma <- matrix(beta[(no.theta*p)+1:(psi*p)], ncol=p)
            g.1 <- g.2 <- term <- third.term <- first.term <- second.term <- NULL
            for(j in 1:p){
                first.term[j] <- -1 * lambda.1 * sqrt(sum((gamma[(((j-1)*psi)+1):(((j-1)*psi)+psi)])^2))
                # first.term[j] <- -1 * lambda.1 * abs(sum(gamma[(((j-1)*psi)+1):(((j-1)*psi)+psi)]))
                second.term[j] <- -1 * lambda.2 * sum((gamma[(((j-1)*psi)+1):(((j-1)*psi)+psi)])^2)
                third.term[j] <- -1 * lambda.3 * sum(abs(theta[,j]))
                term[j] <- first.term[j] + second.term[j] + third.term[j]
                # term[j] <- first.term[j] + second.term[j] - lambda.3 * sum(abs(beta))
            }
            return(sum(term))
        }
        return(log.lik(beta) + log.prior(beta))
    }
    # beta.emp <- c(rep(0, no.theta*p), rep(0, p*psi))
    beta.emp <- c(as.vector(theta.origin), as.vector(gamma.origin))
    beta.map <- optim(beta.emp, fn = log.posterior, #gr = grad.log.posterior, 
                    y.origin = y.origin,
                    # method = "BFGS",
                    method = "CG",
                    # method = "SANN",
                    control = list(fnscale = -1))
    # theta.map <- matrix(beta.map$par[1:(2*p)],nrow=2)
    theta.map <- beta.map$par[1:(no.theta*p)]
    gamma.map <- beta.map$par[(no.theta*p)+1:(psi*p)]

    # df.theta <- data.frame("seq" = seq(1, (no.theta*p)),
    #                 theta.map,
    #                 "theta.true" = as.vector(theta.origin))
    # df.theta$covariate <- factor(as.character(rep(seq(1, 1 + nrow(df.theta) %/% no.theta), each = no.theta, length.out = nrow(df.theta))))
    # df.theta$labels <- factor(1:(no.theta*p))
    # ggplot(df.theta, aes(x = labels)) + ylab("") + 
    # geom_point(aes(y = theta.map, color = covariate), size = 1.5) + 
    # geom_point(aes(y = theta.true, color = "true"), size = 2.5) +
    # labs(title=expression("MAP vs True for"~theta)) + 
    # scale_x_discrete(labels = c("",expression(bold(theta[1])),
    #                             "",expression(bold(theta[2])),
    #                             "",expression(bold(theta[3])),
    #                             "",expression(bold(theta[4])),
    #                             "",expression(bold(theta[5])),
    #                             "",expression(bold(theta[6])),
    #                             "",expression(bold(theta[7])),
    #                             "",expression(bold(theta[8])),
    #                             "",expression(bold(theta[9])),
    #                             "",expression(bold(theta[10]))))+
    # theme(plot.title = element_text(hjust = 0.5, size = 20),
    #         axis.ticks.x = element_blank(),
    #         axis.text.x = element_text(hjust=2),
    #         panel.grid.minor.x = element_blank())
        
    # df <- data.frame("seq" = seq(1, (psi*p)), 
    #                 gamma.map, 
    #                 "gamma.true" = as.vector(gamma.origin))
    # df$covariate <- factor(rep(seq(1, 1 + nrow(df) %/% psi), each = psi, length.out = nrow(df)))
    # df$labels <- factor(1:(psi*p))
    # ggplot(df, aes(x =labels , y = gamma.map, color = covariate)) + 
    #     geom_point() + ylab("") +
    #     # geom_smooth(method="gam") +
    #     geom_point(aes(y = gamma.true, color = "true")) + 
    #     # geom_line(aes(x = seq, y = true, color = "true"), linetype = 2) +
    #     labs(title=expression("MAP vs True for"~gamma)) + 
    #     # ggtitle(expression(atop(paste("MAP vs True for ", bold(gamma))))) +
    #     #   annotate("text", x = seq(0, 330, length.out=10), y = -1, label = beta, colour = "red", size = 10) +
    #     scale_x_discrete(labels = c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", expression(bold(gamma[1])),"", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", 
    #     expression(bold(gamma[2])),"", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", expression(bold(gamma[3])),"", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", expression(bold(gamma[4])),"", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", expression(bold(gamma[5])),"", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", expression(bold(gamma[6])),"", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", expression(bold(gamma[7])),"", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", expression(bold(gamma[8])),"", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", expression(bold(gamma[9])),"", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", expression(bold(gamma[10]))),10)+
    #     theme(plot.title = element_text(hjust = 0.5, size = 20),
    #             axis.ticks.x = element_blank(),
    #             axis.text.x = element_text(hjust=4),
    #             panel.grid.major.x = element_blank())

    # f.nonlinear.new <- f.linear.new <- f.new <- matrix(, nrow = n, ncol=p)
    # for (j in 1:p){
    #     f.linear.new[,j] <- as.matrix(bs.linear[,(((j-1)*no.theta)+1):(((j-1)*no.theta)+no.theta)]) %*% matrix(theta.map, nrow=no.theta)[,j]
    #     f.nonlinear.new[,j] <- bs.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% matrix(gamma.map, ncol = p)[, j]
    #     f.new[1:n, j] <- f.linear.new[,j] + f.nonlinear.new[,j]
    # }
    # new.y <- newalpha <- NULL

    # # set.seed(100)
    # for(i in 1:n){  
    #     temp <- sum(f.new[i,])
    #     newalpha[i] <- exp(temp)
    #     # new.y[i] <- rPareto(1, 1, alpha = newalpha[i])
    # }

    # func.linear.new <- func.nonlinear.new <- func.linear.origin <- func.nonlinear.origin <- func.new <- func.origin <- matrix(, nrow=n, ncol=0)
    # for (j in 1:p){
    #     func.origin <- cbind(func.origin, sort(f.origin[,j]))
    #     func.linear.origin <- cbind(func.linear.origin, sort(f.linear.origin[,j]))
    #     func.nonlinear.origin <- cbind(func.nonlinear.origin, sort(f.nonlinear.origin[,j]))
    #     func.new <- cbind(func.new, sort(f.new[,j]))
    #     func.linear.new <- cbind(func.linear.new, sort(f.linear.new[,j]))
    #     func.nonlinear.new <- cbind(func.nonlinear.new, sort(f.nonlinear.new[,j]))  
    # }
    # # for (j in 1:p){
    # #   func.new <- cbind(func.new, sort(f.new[,j]))
    # # }
    # covariates <- gl(p, n, (p*n))
    # replicate <- gl(2, n, (p*n))
    # func.df <- data.frame(seq = seq(0,1,length.out = n),
    #                         x = as.vector(apply(x.origin, 2, sort, method = "quick")),
    #                         origin=as.vector(func.origin),
    #                         origin.linear=as.vector(func.linear.origin),
    #                         origin.nonlinear=as.vector(func.nonlinear.origin), 
    #                         new=as.vector(func.new),
    #                         new.linear=as.vector(func.linear.new),
    #                         new.nonlinear=as.vector(func.nonlinear.new),
    #                         covariates=covariates, 
    #                         replicate=replicate)

    # ggplot(func.df, aes(x=x, group=interaction(covariates, replicate))) + 
    #     geom_line(aes(y=origin, colour = covariates, linetype = "true")) + 
    #     geom_line(aes(y=new, colour = covariates, linetype = "MAP")) + ylab ("") +
    #     # geom_point(aes(y=origin, shape = replicate)) + geom_point(aes(y=new, shape = replicate)) +
    #     facet_grid(covariates ~ .) + ggtitle("MAP vs True for Smooth Functions") + 
    #     scale_linetype_manual("functions",values=c("MAP"=3,"true"=1)) +
    #     theme(plot.title = element_text(hjust = 0.5, size = 20))

    # ggplot(func.df, aes(x=x, group=interaction(covariates, replicate))) + 
    #     geom_line(aes(y=origin.linear, colour = covariates, linetype = "true")) + 
    #     geom_line(aes(y=new.linear, colour = covariates, linetype = "MAP")) + ylab ("") +
    #     # geom_point(aes(y=origin, shape = replicate)) + geom_point(aes(y=new, shape = replicate)) +
    #     facet_grid(covariates ~ .) + ggtitle("MAP vs True for Linear Components of Smooth Functions") + 
    #     scale_linetype_manual("functions",values=c("MAP"=3,"true"=1)) +
    #     theme(plot.title = element_text(hjust = 0.5, size = 20))

    # ggplot(func.df, aes(x=x, group=interaction(covariates, replicate))) + 
    #     geom_line(aes(y=origin.nonlinear, colour = covariates, linetype = "true")) + 
    #     geom_line(aes(y=new.nonlinear, colour = covariates, linetype = "MAP")) + ylab ("") +
    #     # geom_point(aes(y=origin, shape = replicate)) + geom_point(aes(y=new, shape = replicate)) +
    #     facet_grid(covariates ~ .) + ggtitle("MAP vs True for Nonlinear Components of Smooth Functions") + 
    #     scale_linetype_manual("functions",values=c("MAP"=3,"true"=1)) +
    #     theme(plot.title = element_text(hjust = 0.5, size = 20))


    # data.scenario <- data.frame("x" = c(1:n),
    #                             "constant" = newx,
    #                             "trueAlp" = sort(alp.origin),
    #                             "mapAlp" = sort(newalpha))

    # plt.samp <- ggplot(data = data.scenario, aes(x = constant)) + ylab(expression(alpha(x)))
    # print(plt.samp + 
    #     geom_line(aes(y = trueAlp, col = paste0("True Alpha:",n,"/",psi,"/",threshold)), linewidth = 2.5) + 
    #     geom_line(aes(y = mapAlp, col = paste0("MAP Alpha:",lambda.1,"/",lambda.2,"/",lambda.3)), linewidth = 2.5, linetype = 2) +
    #     labs(col = "") +
    #     scale_color_manual(values = c("#e0b430", "red"))+
    #     theme(text = element_text(size = 15),
    #             legend.position="bottom", legend.key.size = unit(1, 'cm'),
    #             axis.text = element_text(size = 20),
    #             legend.margin=margin(-15,-15,-15,-15),
    #             legend.box.margin=margin(-25,0,20,0)))

    # plot(sort(alp.origin))
    # plot(sort(new.y), sort(y.origin))
    # abline(a=0, b=1, col = "red", lty = 2)
    # rbind(matrix(theta.map, nrow = no.theta, ncol = p), matrix(gamma.map, nrow = psi, ncol = p))

    #### Test Set
    f.nonlinear.origin <- f.linear.origin <- f.origin <- f.nonlinear.new <- f.linear.new <- f.new <- matrix(, nrow = n, ncol=p)
    true.alpha <- new.alpha <- NULL
    for (j in 1:p){
        f.linear.new[,j] <- as.matrix(xholder.linear[,(((j-1)*no.theta)+1):(((j-1)*no.theta)+no.theta)]) %*% matrix(theta.map,nrow=no.theta)[,j]
        f.nonlinear.new[,j] <- xholder.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% matrix(gamma.map, ncol = p)[, j]
        f.new[,j] <- f.linear.new[,j] + f.nonlinear.new[,j]
        f.linear.origin[,j] <- as.matrix(xholder.linear[,(((j-1)*no.theta)+1):(((j-1)*no.theta)+no.theta)]) %*% theta.origin[,j]
        f.nonlinear.origin[,j] <- xholder.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j]
        f.origin[,j] <- f.linear.origin[,j] + f.nonlinear.origin[,j]
    }
        # set.seed(100)
    for(i in 1:n){
        true.alpha[i] <- exp(sum(f.origin[i,]))
        new.alpha[i] <- exp(sum(f.new[i,]))
    }

    func.linear.new <- func.nonlinear.new <- func.linear.origin <- func.nonlinear.origin <- func.new <- func.origin <- matrix(, nrow=n, ncol=0)
    for (j in 1:p){
        func.origin <- cbind(func.origin, f.origin[,j])
        func.linear.origin <- cbind(func.linear.origin, f.linear.origin)
        func.nonlinear.origin <- cbind(func.nonlinear.origin, f.nonlinear.origin)
        func.new <- cbind(func.new, f.new[,j])
        func.linear.new <- cbind(func.linear.new, f.linear.new)
        func.nonlinear.new <- cbind(func.nonlinear.new, f.nonlinear.new)
    }

    covariates <- gl(p, n, (p*n))
    replicate <- gl(2, n, (p*n))
    func.df <- data.frame(seq = seq(0,1,length.out = n), 
                            origin = as.vector(func.origin), 
                            new = as.vector(func.new), 
                            covariates=covariates, 
                            replicate=replicate)

    # ggplot(func.df, aes(x=seq, group=interaction(covariates, replicate))) + 
    #     geom_line(aes(y=origin, colour = covariates, linetype = "true")) + 
    #     geom_line(aes(y=new, colour = covariates, linetype = "MAP")) + ylab ("") +
    #     # geom_point(aes(y=origin, shape = replicate)) + geom_point(aes(y=new, shape = replicate)) +
    #     facet_grid(covariates ~ .) + ggtitle("MAP vs True for Smooth Functions") + 
    #     scale_linetype_manual("functions",values=c("MAP"=3,"true"=1)) +
    #     theme(plot.title = element_text(hjust = 0.5, size = 20))

    data.scenario <- data.frame("x" = c(1:n),
                                "constant" = newx,
                                "trueAlp" = sort(true.alpha),
                                "mapAlp" = sort(new.alpha))

    plt.samp <- ggplot(data = data.scenario, aes(x = constant)) + ylab(expression(alpha(x))) + xlab("")
    print(plt.samp + 
        geom_line(aes(y = trueAlp, col = paste0("True Alpha:",n,"/",psi,"/",threshold)), linewidth = 2.5) + 
        geom_line(aes(y = mapAlp, col = paste0("MAP Alpha:",lambda.1,"/",lambda.3)), linewidth = 2.5, linetype = 2) +
        labs(col = "") +
        theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
        theme(axis.title.x = element_text(size = rel(1.8), angle = 00)) +
        scale_color_manual(values = c("#e0b430", "red"))+
        theme(text = element_text(size = 15),
                legend.position="bottom", legend.key.size = unit(1, 'cm'),
                axis.text = element_text(size = 20),
                legend.margin=margin(-15,-15,-15,-15),
                legend.box.margin=margin(-25,0,20,0)))

    # Randomized quantile residuals
    # r <- matrix(NA, nrow = n, ncol = 20)
    # # beta <- as.matrix(mcmc[[1]])[, 1:7] 
    # T <- 20
    # for(i in 1:n){
    #     for(t in 1:T){
    #         r[i, t] <- qnorm(pPareto(y.origin[i], 1, alpha = newalpha[i]))
    #     }
    # }
    # lgrid <- n
    # grid <- qnorm(ppoints(lgrid))
    # # qqnorm(r[, 1])
    # # points(grid, quantile(r[, 1], ppoints(lgrid), type = 2), 
    # #     xlim = c(-3, 3), col = "red")
    # traj <- matrix(NA, nrow = T, ncol = lgrid)
    # for (t in 1:T){
    # traj[t, ] <- quantile(r[, t], ppoints(lgrid), type = 2)
    # }
    # l.band <- apply(traj, 2, quantile, prob = 0.025)
    # trajhat <- apply(traj, 2, quantile, prob = 0.5)
    # u.band <- apply(traj, 2, quantile, prob = 0.975)

    # ggplot(data = data.frame(grid = grid, l.band = l.band, trajhat = trajhat, 
    #                         u.band = u.band)) + 
    # #geom_ribbon(aes(x = grid, ymin = l.band, ymax = u.band), 
    # #           color = "lightgrey", fill = "lightgrey",
    # #          alpha = 0.4, linetype = "dashed") + 
    # geom_ribbon(aes(x = grid, ymin = l.band, ymax = u.band), 
    #             color = "lightgrey", fill = "lightgrey",
    #             alpha = 0.4, linetype = "dashed") + 
    # geom_line(aes(x = grid, y = trajhat), linetype = "dashed", linewidth = 1.2) + 
    # geom_abline(intercept = 0, slope = 1, linewidth = 1.2) + 
    # labs(x = "Theoretical quantiles", y = "Sample quantiles") + 
    # # theme(text = element_text(size = 30)) + 
    # coord_fixed(xlim = c(-3, 3),  
    #             ylim = c(-3, 3)) 
    # ggsave("../../Figures/residuals_application.pdf")
    # ggsave("../Results/residuals_application_IWSM.pdf", width=6, height=6)


#--------------------------------------------------------------------------------------

  alpha.container[,mm] <- data.scenario$mapAlp
  gamma.container[,mm] <- gamma.map
#   theta0.container[,mm] <- matrix(theta.map, nrow=2)[1,]
  theta.container[,mm] <- theta.map
  linear.container[[mm]] <- func.linear.new
  nonlinear.container[[mm]] <- func.nonlinear.new
  f.container[[mm]] <- func.new
#   f.container[,mm] <- as.vector(func.new)
}
alpha.container$x <- seq(0,1, length.out = n)
alpha.container$true <- data.scenario$trueAlp
alpha.container <- cbind(alpha.container, t(apply(alpha.container[,1:simul.no], 1, quantile, c(0.05, .5, .95))))
colnames(alpha.container)[(dim(alpha.container)[2]-2):(dim(alpha.container)[2])] <- c("q1","q2","q3")
alpha.container$mean <- rowMeans(alpha.container[,1:simul.no])
alpha.container <- as.data.frame(alpha.container)
plt <- ggplot(data = alpha.container, aes(x = x)) + ylab(expression(alpha(x))) + xlab("")
for(i in 1:simul.no){
  # print(.data[[names(data.scenario)[i]]])
  plt <- plt + geom_line(aes(y = .data[[names(alpha.container)[i]]]), alpha = 0.4,linewidth = 0.7)
  # plt <- plt + geom_line(aes(y = .data[[names(data.scenario)[i]]]))
}
print(plt + geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.5) + 
        geom_line(aes(y=true, col = "True"), linewidth = 1.8) + 
        # geom_line(aes(y = post.check, col=paste0("Simulated Alpha: ",n,"/",psi,"/",threshold)), linewidth = 1.5) +
        geom_line(aes(y=mean, col = "Mean"), linewidth = 1.8, linetype = 2) +
        theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
        theme(axis.title.x = element_text(size = rel(1.8), angle = 00)) +
        labs(col = "") + ylim(0, 150) + 
        scale_color_manual(values = c("#e0b430", "red"))+
        theme(text = element_text(size = 15),
                legend.position="bottom", legend.key.size = unit(1, 'cm'),
                axis.text = element_text(size = 20),
                legend.margin=margin(-15,-15,-15,-15),
                legend.box.margin=margin(-25,0,20,0)))

ggsave(paste0("../Laboratory/Simulation/BayesianFusedLasso/mcresults/map_alpha_n", simul.no,".pdf"), width=10)

plt <- ggplot(data = alpha.container[900:1000,], aes(x = x)) + ylab(expression(alpha(x))) + xlab("")
for(i in 1:simul.no){
  # print(.data[[names(data.scenario)[i]]])
  plt <- plt + geom_line(aes(y = .data[[names(alpha.container)[i]]]), alpha = 0.4,linewidth = 0.7)
  # plt <- plt + geom_line(aes(y = .data[[names(data.scenario)[i]]]))
}
print(plt + geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.5) + 
        geom_line(aes(y=true, col = "True"), linewidth = 1.8) + 
        geom_line(aes(y=mean, col = "Mean"), linewidth = 1.8, linetype = 2) +
        # geom_line(aes(y = post.check, col=paste0("Simulated Alpha: ",n,"/",psi,"/",threshold)), linewidth = 1.5) +
        theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
        theme(axis.title.x = element_text(size = rel(1.8), angle = 00)) +
        labs(col = "") + ylim(0, 150) +
        scale_color_manual(values = c("#e0b430", "red"))+
        theme(text = element_text(size = 15),
                legend.position="bottom", legend.key.size = unit(1, 'cm'),
                axis.text = element_text(size = 20),
                legend.margin=margin(-15,-15,-15,-15),
                legend.box.margin=margin(-25,0,20,0)))
    
ggsave(paste0("../Laboratory/Simulation/BayesianFusedLasso/mcresults/map_alpha_90q_n", simul.no,".pdf"), width=10)

resg <- gather(theta.container,
               key = "group",
               names(theta.container),
               value = "values")
resg$group1 <- factor(rep(1:(no.theta*p), simul.no))
resg$group2 <- factor(rep(1:p, each = no.theta))
# resg$group2 <- factor(rep(c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10"), each = no.theta))
# factor(c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10"))
# resg$group2 <- factor(rep(1:(no.theta*p), each = (no.theta*p), length.out = nrow(no.theta*p)))
somelines <- data.frame(value=c(as.vector(theta.origin)),boxplot.nr=c(1:(no.theta*p)))
ggplot(resg, aes(group=group1, x = group1, y = values, fill=group2)) + ylim(-1,1) + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) +
  geom_boxplot() + #coord_cartesian(ylim=c(-1,1))+
  geom_segment(data=somelines,aes(x=boxplot.nr-0.5,xend=boxplot.nr+0.5, 
                                  y=value,yend=value),inherit.aes=FALSE,color="red",linewidth=1.5)+
  # facet_wrap( ~ group2, labeller = label_parsed) +
  scale_x_discrete(labels = c("", expression(theta[1]), 
                              "", expression(theta[2]),
                              "", expression(theta[3]),
                              "", expression(theta[4]),
                              "", expression(theta[5]),
                              "", expression(theta[6]),
                              "", expression(theta[7]),
                              "", expression(theta[8]),
                              "", expression(theta[9]),
                              "", expression(theta[10]))) +
  labs(x = "", y = "") + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_text(hjust=0.8),
        axis.text = element_text(size = 30),
        legend.title = element_blank(),
        legend.text = element_text(size=30))

ggsave(paste0("../Laboratory/Simulation/BayesianFusedLasso/mcresults/map_theta_n", simul.no,".pdf"), width=13)

resg <- gather(gamma.container,
               key = "group",
               names(gamma.container),
               value = "values")
resg$group1 <- factor(rep(1:(psi*p), simul.no))
resg$group2 <- factor(rep(1:p, each = psi))
# resg$group2 <- factor(rep(c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10"), each = psi))

# resg$group2 <- factor(rep(1:(no.theta*p), each = (no.theta*p), length.out = nrow(no.theta*p)))
somelines <- data.frame(value=c(as.vector(gamma.origin)),
                        boxplot.nr=c(1:(psi*p)),
                        covariate = factor(rep(1:p, each= psi)))
ggplot(resg, aes(group=group1, x = group1, y = values, fill=group2)) + ylim(-1.1,1.1) + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) +
  geom_boxplot() + #coord_cartesian(ylim=c(-1,1))+
  geom_segment(data=somelines,aes(x=boxplot.nr-0.5,xend=boxplot.nr+0.5, 
                                  y=value,yend=value),inherit.aes=FALSE,
                                  color="red",
                                  linewidth=1.5)+
  # facet_wrap( ~ group2, labeller = label_parsed) +
  scale_x_discrete(labels = c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", expression(bold(gamma[1])),
                        "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", expression(bold(gamma[2])), "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", expression(bold(gamma[3])),"", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", expression(bold(gamma[4])),"", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", expression(bold(gamma[5])),"", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", expression(bold(gamma[6])),"", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", expression(bold(gamma[7])),
                        "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", expression(bold(gamma[8])),
                        "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", expression(bold(gamma[9])),
                        "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", expression(bold(gamma[10])),
                        "", "", "", "", "", "", "", "", "", "", "", "", "")) +
  labs(x = "", y = "") + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_text(hjust=1.75),
        axis.text = element_text(size = 30),
        legend.title = element_blank(),
        legend.text = element_text(size=30),        
        panel.grid.major.x = element_blank())

ggsave(paste0("../Laboratory/Simulation/BayesianFusedLasso/mcresults/map_gamma_n", simul.no,".pdf"), width=13)

# bb <- as.data.frame(f.container)
# head(bb)
quant <- apply(simplify2array(f.container), 1:2, quantile, c(0.05, 0.5, 0.95))
q1 <- quant[1, , ]
q2 <- quant[2, , ]
q3 <- quant[3, , ]
names(f.container) <- paste("s", 1:simul.no, sep = "")
f.container <- c(f.container, list(q1=q1, q2 = q2, q3 =q3, true = func.origin, x=xholder, colour = matrix(rep(rainbow(p), each=n),nrow=n)))
str(f.container)
# dresults <- lapply(f.container, as.data.frame) %>% bind_rows(.id = "id")
dresults <- bind_rows(f.container, .id = "data_frame")
plts <- list()
for (i in 1:p){
    if(i != p){
        plt <- ggplot(dresults, aes(x=x[,!! i], y=true[, !! i])) + 
            geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 1.5) +
            geom_ribbon(aes(ymin = q1[, !! i], ymax = q3[, !! i]), alpha = 0.5) + xlab("") + ylab("") + ylim(min(dresults$q1[,i]) -0.1, max(dresults$q3[,i])+0.1) +
            geom_line(aes(col=colour[,!! i]), linewidth = 1.8) +
            geom_line(aes(y=q2[, !! i], col = "Mean"), linewidth = 1.8, linetype = 2) + 
            scale_colour_manual(values = c(as.character(hue_pal()(p)[i]), as.character(hue_pal()(10)[i]))) + 
            theme(#axis.ticks = element_blank(),
                legend.position = "none",
                # legend.title = element_blank(),
                # legend.text = element_text(size=33),
                # strip.text = element_blank(),
                axis.text = element_blank(),
                plot.margin = unit(c(0, 3.5, -8, 3.5), "pt"))
    }
    else{
        plt <- ggplot(dresults, aes(x=x[,!! i], y=true[, !! i])) + 
            geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 1.5) +
            geom_ribbon(aes(ymin = q1[, !! i], ymax = q3[, !! i]), alpha = 0.5) + xlab("Smooth Functions") + ylab("") +  ylim(min(dresults$q1[,i]) -0.1, max(dresults$q3[,i])+0.1) +
            geom_line(aes(col=colour[,!! i]), linewidth = 1.8) +
            geom_line(aes(y=q2[, !! i], col = "Mean"), linewidth = 1.8, linetype = 2) + 
            scale_colour_manual(values = c(as.character(hue_pal()(p)[i]), as.character(hue_pal()(p)[i]))) + 
            theme(legend.position = "none",
                    plot.margin = unit(c(1, 3.5, 0, 3.5), "pt"),
                    axis.text.y = element_blank())#,
                    # axis.title.x = element_text(size = 33))
    }
    plts[[i]] <- plt
}
# plts
# grid.arrange(grobs = plts, nrow = p)
g <- grid.arrange(grobs = plts, nrow = p)
g# g <-arrangeGrob(plts, nrow = p)
ggsave(file = paste0("../Laboratory/Simulation/BayesianFusedLasso/mcresults/map_smooth_", simul.no,".pdf"), g, width=10, height = 16)
