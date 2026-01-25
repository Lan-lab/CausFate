# This file includes code from:

# Author: Minzhe Zhang
# Title: geneck
# Source: https://github.com/Minzhe/geneck
# License: Creative Commons Attribution 3.0 Unported (CC BY 3.0)
#          https://creativecommons.org/licenses/by/3.0/

network.cmi2ni <- function(expr.data, lambda) {
  if (lambda <= 0) {
    stop("Input error: parameter lambda for cmi2ni should be larger than 0.")
  }
  p <- ncol(expr.data)
  n <- nrow(expr.data)
  gene.index <- colnames(expr.data)

  expr.mat <- scale(as.matrix(expr.data), center = TRUE, scale = FALSE)

  out <- cmi2ni(t(expr.mat), lambda)

  est_edge <- which(out$G == 1, T)
  est_edge <- est_edge[est_edge[, 1] < est_edge[, 2], ]
  if (length(est_edge) == 2) {
    est_edge <- matrix(est_edge, 1, 2)
  }

  est_edge <- as.data.frame(est_edge)
  if (nrow(est_edge) == 0) {
    (return(data.frame(
      node1 = character(),
      node2 = character(),
      strength = double()
    )))
  }
  colnames(est_edge) <- c("node1", "node2")
  est_edge$strength <- 0
  for (i in 1:nrow(est_edge)) {
    est_edge$strength[i] <- out$Gval[est_edge$node1[i], est_edge$node2[i]]
  }

  est_edge$strength <- signif(est_edge$strength, 4)
  est_edge <- est_edge[with(est_edge, order(-abs(strength))), ]
  est_edge[, 1] <- gene.index[est_edge[, 1]]
  est_edge[, 2] <- gene.index[est_edge[, 2]]

  return(est_edge)
}

cmi <- function(v1, v2, nargin = 2, vcs = NULL) {
  if (nargin == 2) {
    c1 <- stats::var(v1)

    c2 <- stats::var(v2)

    c3 <- det(stats::cov(t(rbind(v1, v2))))

    cmiv <- 0.5 * log(c1 * c2 / c3)
  } else if (nargin == 3) {
    c1 <- det(stats::cov(t(rbind(v1, vcs))))

    c2 <- det(stats::cov(t(rbind(v2, vcs))))

    if (is.null(dim(vcs))) {
      c3 <- stats::var(vcs)
    } else {
      c3 <- det(stats::cov(t(vcs)))
    }
    c4 <- det(stats::cov(t(rbind(v1, v2, vcs))))

    cmiv <- 0.5 * log((c1 * c2) / (c3 * c4))
  }

  # if(nargin==2)
  # cat('\n values: ',c1,c2,c3,cmiv,sep=' ')
  # if(nargin==3)
  # cat('\n values: ',c1,c2,c3,c4,cmiv,sep=' ')

  if (cmiv == Inf) {
    cmiv <- 1e10
  }

  return(cmiv)
}

cmi2ni <- function(dat, lambda, order0 = 0, nargin = 2) {
  # function [G,Gval,order]=pca_cmi(data,lamda,order0)

  n_gene <- nrow(dat)
  G <- matrix(1, n_gene, n_gene)
  G <- t(as.matrix(Matrix::tril(G, -1)))

  G <- G + t(G)

  Gval <- G

  order <- -1

  t <- 0

  while (t == 0) {
    order <- order + 1

    if (nargin == 3) {
      if (order > order0) {
        G <- t(as.matrix(Matrix::tril(G, -1)))

        Gval <- t(as.matrix(Matrix::tril(G, -1)))

        order <- order - 1
        # The value of order is the last order of pc algorith
        return(list(
          G = G,
          Gval = Gval,
          order = order
        ))
      }
    }

    res <- edgereduce(G, Gval, order, dat, t, lambda)

    G <- res$G
    Gval <- res$Gval
    t <- res$t

    if (t == 0) {
      # cat('\n No edge is reduce! Algorithm  finished!')

      break
    } else {
      t <- 0
    }
  }

  G <- t(as.matrix(Matrix::tril(G, -1)))

  Gval <- t(as.matrix(Matrix::tril(Gval, -1)))

  order <- order - 1
  # The value of order is the last order of pc algorith
  return(list(
    G = G,
    Gval = Gval,
    order = order
  ))
}

cas <- function(x, y, z) {
  # x=rand(10,1)';y=rand(10,1)';z=rand(10,2)';
  if (is.null(dim(z))) {
    n1 <- 1
  } else {
    n1 <- nrow(z)
  }
  n <- n1 + 2

  Cov <- as.matrix(stats::var(x))

  Covm <- stats::cov(t(rbind(x, y, z)))

  Covm1 <- stats::cov(t(rbind(x, z)))

  InvCov <- solve(Cov)

  InvCovm <- solve(Covm)

  InvCovm1 <- solve(Covm1)

  C11 <- InvCovm1[1, 1]

  C12 <- 0

  C13 <- InvCovm1[1, 2:(1 + n1)]

  C23 <- InvCovm[2, 3:(2 + n1)] -
    InvCovm[1, 2] *
      (1 / (InvCovm[1, 1] - InvCovm1[1, 1] + InvCov[1, 1])) *
      (InvCovm[1, 3:(2 + n1)] - InvCovm1[1, 2:(1 + n1)])

  C22 <- InvCovm[2, 2] -
    InvCovm[1, 2]^2 * (1 / (InvCovm[1, 1] - InvCovm1[1, 1] + InvCov[1, 1]))

  C33 <- InvCovm[3:(2 + n1), 3:(2 + n1)] -
    (1 / (InvCovm[1, 1] - InvCovm1[1, 1] + InvCov[1, 1])) *
      ((InvCovm[1, 3:(2 + n1)] - InvCovm1[1, 2:(1 + n1)]) %o%
        (InvCovm[
          1,
          3:(2 +
            n1)
        ] -
          InvCovm1[1, 2:(1 + n1)]))

  InvC <- rbind(c(C11, C12, C13), c(C12, C22, C23), cbind(C13, C23, C33))

  # %C = inv(InvC);

  C0 <- Cov[1, 1] * (InvCovm[1, 1] - InvCovm1[1, 1] + InvCov[1, 1])

  CS <- 0.5 * (sum(diag(InvC %*% Covm)) + log(C0) - n)

  return(CS)
}

edgereduce <- function(G, Gval, order, dat, t, lambda, nargin = 2) {
  if (order == 0) {
    for (i in 1:(nrow(G) - 1)) {
      for (j in (i + 1):nrow(G)) {
        if (G[i, j] != 0) {
          cmiv <- cmi(dat[i, ], dat[j, ])

          Gval[i, j] <- cmiv

          Gval[j, i] <- cmiv

          if (cmiv < lambda) {
            G[i, j] <- 0

            G[j, i] <- 0
          }
        }
      }
    }
    t <- t + 1
  } else {
    for (i in 1:(nrow(G) - 1)) {
      for (j in (i + 1):nrow(G)) {
        if (G[i, j] != 0) {
          adj <- numeric()

          for (k in 1:nrow(G)) {
            if (G[i, k] != 0 && G[j, k] != 0) {
              adj <- c(adj, k)
            }
          }
          if (length(adj) >= order) {
            if (length(adj) == 1) {
              combntnslist <- adj
              combntnsrow <- 1
            } else {
              combntnslist <- t(utils::combn(adj, order))

              combntnsrow <- nrow(combntnslist)
            }
            cmiv <- 0

            v1 <- dat[i, ]

            v2 <- dat[j, ]

            if (combntnsrow == 1) {
              vcs <- dat[combntnslist, ]
              # cat('\n ',adj,'adj',combntnsrow,k,i,j,' list ', combntnslist,' dim ', nrow(vcs), ncol(vcs),sep=' ')
              cmiv <- MI2(v1, v2, vcs)
            } else {
              for (k in 1:combntnsrow) {
                vcs <- dat[combntnslist[k, ], ]

                # cat('\n ',adj,'adj',combntnsrow,k,i,j,' list ', combntnslist[k,],' dim ', nrow(vcs), ncol(vcs),sep=' ')
                a <- MI2(v1, v2, vcs)

                cmiv <- max(cmiv, a)
              }
            }
            Gval[i, j] <- cmiv

            Gval[j, i] <- cmiv

            if (cmiv < lambda) {
              G[i, j] <- 0

              G[j, i] <- 0
            }
            t <- t + 1
          }
        }
      }
    }
  }
  return(list(G = G, Gval = Gval, t = t))
}

MI2 <- function(x, y, z) {
  r_dmi <- (cas(x, y, z) + cas(y, x, z)) / 2

  return(r_dmi)
}

###                        space                           ###
### ====================================================== ###
# The following R scripts are functions to use space to construct gene network.

network.space <- function(expr.data, alpha) {
  if (alpha <= 0) {
    stop("Input error: parameter alpha for space should be larger than 0.")
  }
  p <- ncol(expr.data)
  n <- nrow(expr.data)
  gene.index <- colnames(expr.data)

  expr.mat <- scale(as.matrix(expr.data), center = TRUE, scale = FALSE)

  out <- space::space.joint(expr.mat, lam1 = alpha * n, iter = 5)
  est_edge <- which(abs(out$ParCor) > 0, T)
  est_edge <- est_edge[est_edge[, 1] < est_edge[, 2], ]
  if (length(est_edge) == 2) {
    est_edge <- matrix(est_edge, 1, 2)
  }

  est_edge <- as.data.frame(est_edge)
  colnames(est_edge) <- c("node1", "node2")
  est_edge$p.corr <- 0
  for (i in 1:nrow(est_edge)) {
    est_edge$p.corr[i] <- out$ParCor[est_edge$node1[i], est_edge$node2[i]]
  }

  est_edge$p.corr <- signif(est_edge$p.corr, 4)
  est_edge <- est_edge[with(est_edge, order(-abs(p.corr))), ]
  est_edge[, 1] <- gene.index[est_edge[, 1]]
  est_edge[, 2] <- gene.index[est_edge[, 2]]

  return(est_edge)

  bic_val <- sum(log(n / out$sig.fit)) + log(n) / n * nrow(est_edge) * 2
}
