inla.barrier.fem.copy <-
function (mesh, barrier.triangles, Omega = NULL)
{
    stopifnot(inherits(mesh, "inla.mesh"))
    if (missing(barrier.triangles) && is.null(Omega))
        stop("Input barrier triangles")
    if (missing(barrier.triangles)) {
    }
    else {
        barrier.triangles <- unique(barrier.triangles)
        t <- length(mesh$graph$tv[, 1])
        remaining <- setdiff(1:t, barrier.triangles)
        if (!is.null(Omega))
            warning("Omega is replaced by barrier.triangles")
        Omega <- list(remaining, barrier.triangles)
    }
    dt.fem.white <- function(mesh, subdomain) {
        Ck <- rep(0, mesh$n)
        for (t in subdomain) {
            px <- mesh$graph$tv[t, ]
            temp <- mesh$loc[px, ]
            p1 <- t(t(temp[1, c(1, 2)]))
            p2 <- t(t(temp[2, c(1, 2)]))
            p3 <- t(t(temp[3, c(1, 2)]))
            Ts <- cbind(p2 - p1, p3 - p1)
            area <- abs(det(Ts)) * 0.5
            for (i in 1:3) {
                Ck[px[i]] <- Ck[px[i]] + area
            }
        }
        return(Ck)
    }
    dt.fem.identity <- function(mesh) {
        len <- length(mesh$graph$tv[, 1])
        index.i <- rep(0, len * 6)
        index.j <- rep(0, len * 6)
        Aij <- rep(0, len * 6)
        counter <- 1
        for (t in 1:len) {
            px <- mesh$graph$tv[t, ]
            temp <- mesh$loc[px, ]
            p1 <- t(t(temp[1, c(1, 2)]))
            p2 <- t(t(temp[2, c(1, 2)]))
            p3 <- t(t(temp[3, c(1, 2)]))
            Ts <- cbind(p2 - p1, p3 - p1)
            twiceArea <- abs(det(Ts))
            for (i in 1:3) {
                index.i[counter] <- px[i]
                index.j[counter] <- px[i]
                Aij[counter] <- (twiceArea) * 1/12
                counter <- counter + 1
            }
            for (i in 1:2) {
                for (j in (i + 1):3) {
                  index.i[counter] <- px[i]
                  index.j[counter] <- px[j]
                  Aij[counter] <- (twiceArea) * 1/24
                  counter <- counter + 1
                  index.i[counter] <- px[j]
                  index.j[counter] <- px[i]
                  Aij[counter] <- (twiceArea) * 1/24
                  counter <- counter + 1
                }
            }
        }
        I <- Matrix::sparseMatrix(i = index.i, j = index.j, x = Aij,
            dims = c(mesh$n, mesh$n), repr = "T")
        return(I)
    }
    dt.fem.laplace <- function(mesh, subdomain) {
        Nphix <- rbind(c(-1, -1), c(1, 0), c(0, 1))
        len <- length(subdomain)
        index.i <- rep(0, len * 9)
        index.j <- rep(0, len * 9)
        Aij <- rep(0, len * 9)
        counter <- 1
        for (tri in subdomain) {
            px <- mesh$graph$tv[tri, ]
            temp <- mesh$loc[px, ]
            p1 <- t(t(temp[1, c(1, 2)]))
            p2 <- t(t(temp[2, c(1, 2)]))
            p3 <- t(t(temp[3, c(1, 2)]))
            Ts <- cbind(p2 - p1, p3 - p1)
            TTTinv <- solve(t(Ts) %*% Ts)
            area <- abs(det(Ts)) * 0.5
            for (k in 1:3) {
                for (m in 1:3) {
                  tmp <- (3 * m + k - 4) * length(subdomain)
                  index.i[(tmp + counter)] <- px[k]
                  index.j[(tmp + counter)] <- px[m]
                  Aij[(tmp + counter)] <- area * Nphix[k, c(1,
                    2)] %*% TTTinv %*% as.matrix(Nphix[m, c(1,
                    2)])
                }
            }
            counter <- counter + 1
        }
        Dk <- Matrix::sparseMatrix(i = index.i, j = index.j, x = Aij,
            dims = c(mesh$n, mesh$n), repr = "T")
        return(Dk)
    }
    xi <- length(Omega)
    #if (require(INLAspacetime)) {
    #    warning("Using implementation from the `INLAspacetime` package")
    #    fem <- INLAspacetime::mesh2fem.barrier(mesh = mesh, barrier.triangles = Omega[[2L]])
    #}
    #else {
        #warning(paste("Please install the `INLAspacetime` package\n",
        #    "which contains an implementation that runs faster!"))
        fem <- list()
        fem$I <- dt.fem.identity(mesh)
        fem$D <- list()
        fem$C <- list()
        for (k in 1:xi) {
            fem$D[[k]] <- dt.fem.laplace(mesh, Omega[[k]])
        }
        for (k in 1:xi) {
            fem$C[[k]] <- dt.fem.white(mesh, Omega[[k]])
        }
        fem$hdim <- xi
    #}
    return(fem)
}
