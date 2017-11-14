# ------------------------------------------------------------------------------
# mnnd function after Mean Nearest Neighbor Distances
# ------------------------------------------------------------------------------
#' Calculates the Mean Nearest Neighbor Distance
#'
#' Function mnnd accepts a xy-coordinate (xy1) or a pair of
#' xy-coordinates (xy1 and xy2). The arguments xy1 and xy2
#' should be a matrix or a data-frame R objects consist of x and y values.
#' Function mnnd then returns the mean nearest neighbor distance based on
#' Clark and Evans (1954) Ecology: 35(4): 445-453.
#'
#' @importFrom stats dist
#' @param xy1 A set of coordinates
#' @param xy2 A set of coordinates
#' @return The Mean Nearest Neighbor Distance
#' @examples
#' coords1 <- matrix(sample(1:10,10),nrow=5,ncol=2)
#' cooords2 <- matrix(sample(1:10,10),nrow=5,ncol=2)
#' mnnd(coords1,cooords2)
#' @export
mnnd <- function(xy1, xy2 = NULL) {
    nd <- c()
    if (is.null(xy2)) {
        d <- as.matrix(dist(xy1))
        for (i in 1:nrow(xy1)) nd <- rbind(nd, min(d[i, -i]))
    } else {
        for (i in 1:nrow(xy1)) nd <- c(nd, min(as.matrix(dist(rbind(xy1[i, ],
                                                                xy2)))[1, -1]))
    }
    mean(nd)
}


# ------------------------------------------------------------------------------
# smnnd function after Standardized Mean Nearest Neighbor Distances
# ------------------------------------------------------------------------------
#' Calculates the Mean Nearest Neighbor Distances
#'
#' Function smnnd accepts a xy-coordinate (xy1) or a pair of xy-coordinates
#' (xy1 and xy2) and the area of the arena in which those spatial points are
#' distributed.  The arguments xy1 and xy2 should be a matrix or a data-frame R
#' objects consist of x and y values. The unit of the area should be the same as
#' that for x and y values.
#' Function smnnd then returns the mean nearest neighbor distance
#' standardized with the Poisson expectation based on Clark and Evans (1954)
#' Ecology: 35(4): 445-453: E(nnd) = 1/sqrt(n)/2, where n is the number of
#' points within a square area with the sides range [0,1].  This expectation,
#' however, underestimates the distances when the number of points are very low
#' (<10), so a tentative correction term: 1/n/4, was added.  Standardization is
#' achieved by subtracting the corrected E(nnd) from the observed mean nearest
#' neighbor distance rescaled with sqrt(area).
#'
#' @param xy1 A set of coordinates
#' @param xy2 A set of coordinates
#' @param area An area in the same units as xy
#' @return The Standardized Mean Nearest Neighbor Distance
#' @examples
#' coords1 <- matrix(sample(1:10,10),nrow=5,ncol=2)
#' coords2 <- matrix(sample(1:10,10),nrow=5,ncol=2)
#' smnnd(coords1,coords2,100)
#' smnnd(coords2,coords1,100)
#' @export
smnnd <- function(xy1, xy2 = NULL, area) {
    if (is.null(xy2)) {
        mnnd(xy1) / sqrt(area) - 1 / sqrt(nrow(xy1)) / 2 - 1 / nrow(xy1) / 4

    } else {
        mnnd(xy1, xy2) / sqrt(area) - 1 / sqrt(nrow(xy2) + 1) / 2 - 1 /
        (nrow(xy2) + 1) / 4
    }
}


# ------------------------------------------------------------------------------
# csp function after Colony Site Persistence
# ------------------------------------------------------------------------------
#' Calculates the Colony Site Persistence Index
#'
#' Function csf accepts a pair of xy-coordinates, representing colony
#' locations, at two different times (xy1 and xy2) and the area of the in which
#' those colony locations are distributed.  The arguments xy1 and xy2 should be
#' a matrix or a data-frame R objects consist of x and y values.  The unit of
#' the area should be the same as that for x and y values.  Function csf
#' then returns the colony site persistence index of two sets of colonies for two
#' consecutive points in time based on Carrasco and Toquenaga (in preparation).
#' The csf index averages the Standardized Mean Nearest Neighbour Distance
#' (smnnd) between the smnnd calculated from the xy1 colonies to the xy2
#' colonies, and the smnnd calculated from the xy2 colonies to the xy1 colonies.
#' A negative sign is added so that increasing values of csfincex represent
#' increases in colony site persistence.
#'
#' @param xy1 A set of coordinates at time 1
#' @param xy2 A set of coordinates at time 2
#' @param area An area in the same units as xy
#' @return The Colony Site Persistence Index
#' @examples
#' colonies1 <- matrix(sample(1:10,10),nrow=5,ncol=2)
#' colonies2 <- matrix(sample(1:10,10),nrow=5,ncol=2)
#' csf(colonies1,colonies2,100)
#' @export
csp <- function(xy1, xy2, area) {
    -(smnnd(xy1, xy2, area) + smnnd(xy2, xy1, area)) / 2
}


# ------------------------------------------------------------------------------
# cspscaled function after Colony Site Persistence
# ------------------------------------------------------------------------------
#' Calculates the Scaled Colony Site Persistence Index
#'
#' Function csf accepts a pair of xy-coordinates, representing colony
#' locations, at two different times (xy1 and xy2) and the area of the in which
#' those colony locations are distributed.  The arguments xy1 and xy2 should be
#' a matrix or a data-frame R objects consist of x and y values.  The unit of
#' the area should be the same as that for x and y values.  Function csf
#' then returns the colony site persistence index of two sets of colonies for two
#' consecutive points in time based on Carrasco and Toquenaga (in preparation).
#' The csf index averages the Standardized Mean Nearest Neighbour Distance
#' (smnnd) between the smnnd calculated from the xy1 colonies to the xy2
#' colonies, and the smnnd calculated from the xy2 colonies to the xy1 colonies.
#' A negative sign is added so that increasing values of csfincex represent
#' increases in colony site persistence.
#'
#' @param xy1 A set of coordinates at time 1
#' @param xy2 A set of coordinates at time 2
#' @param area An area in the same units as xy
#' @return The Colony Site Persistence Index
#' @examples
#' colonies1 <- matrix(sample(1:10,10),nrow=5,ncol=2)
#' colonies2 <- matrix(sample(1:10,10),nrow=5,ncol=2)
#' csf(colonies1,colonies2,100)
#' @export
cspscaled <- function(xy1, xy2, area) {
  -(smnnd(xy1, xy2, area) + smnnd(xy2, xy1, area)) / 2
}
