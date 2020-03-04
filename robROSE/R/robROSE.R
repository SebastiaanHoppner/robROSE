robROSE <- function (formula, data, r, dup_size, alpha = 0.5, const = 1, seed = NULL) {
  # match call
  call <- match.call()

  # check inputs
  if (missing(formula)) {
    stop("argument 'formula' is missing, with no default")
  }
  if (missing(data)) {
    stop("argument 'data' is missing, with no default")
  }
  if (missing(r)) {
    r <- NULL
  }
  if (missing(dup_size)) {
    dup_size <- NULL
  }
  if (is.null(missing) & is.null(dup_size)) {
    stop("both arguments 'r' and 'dup_size' are missing, with no default. One has to be specified.")
  }
  if (!is.null(r) & !is.null(dup_size)) {
    stop("either 'r' or 'dup_size' must be specfied, not both")
  }
  dup_size <- as.integer(dup_size)

  # adjust formula if needed
  formula_original <- formula
  formula <- terms(formula, data = data)
  vars <- attr(formula, "variables")
  vars <- sapply(vars, function (x) paste(deparse(x, width.cutoff = 500), collapse = ' '))[-1L]
  vars <- sub("*.*[(/]", "", vars) # remove all characters before either "(" or "/"
  vars <- sub("['^')].*", "", vars) # remove all characters after either "^" or ")"
  vars <- unique(vars)
  formula <- as.formula(paste(vars[1], "~", paste(vars[-1], collapse = "+")))
  attr(formula, "variables") <- vars

  if (formula_original[[3]] != "." & eval(formula) != formula_original) {
    stop("Transformations of variables are not allowed.")
  }

  # extract response variable and model matrix
  mf <- model.frame(formula, data)
  colNames <- rownames(attributes(attributes(mf)$terms)$factors)
  y <- mf[, 1]
  X <- mf[, -1]

  n <- NROW(X)
  d <- NCOL(X)

  classy <- class(y)
  y <- factor(y)
  Taby <- table(y)
  classX <- sapply(as.data.frame(X), class)

  # check input data
  if (n < 2) {
    stop("Too few observations.")
  }
  if (length(Taby) > 2) {
    stop("The response variable must have 2 levels.")
  } else if (length(Taby) == 1) {
    stop("The response variable has only one class.")
  }
  if (any(is.na(pmatch(classX, c("numeric", "integer", "factor"), duplicates.ok = TRUE)))) {
    stop("The current implementation of robROSE handles only continuous and categorical variables.")
  }

  # select minority and majority samples
  majoY <- levels(y)[which.max(Taby)]
  minoY <- levels(y)[which.min(Taby)]
  ind.mayo <- which(y == majoY)
  ind.mino <- which(y == minoY)

  # select numeric variables
  id.num <- which(classX == "numeric" | classX == "integer")
  d.num <- d - length(which(classX == "factor"))

  # apply MCD on numeric data of minority samples
  if (!is.null(seed)) set.seed(seed)
  covest <- covMcd(X[ind.mino, id.num], alpha = alpha, nsamp = 500)
  names(covest$mah) <- ind.mino

  # find "clean" and "outlying" (i.e. outcast) minority samples
  threshold <- qchisq(0.999, df = d.num)
  ind.mino.clean <- ind.mino[which(covest$mah <  threshold)]
  ind.mino.out   <- ind.mino[which(covest$mah >= threshold)]
  n.mino.clean <- length(ind.mino.clean)

  # number of new minority samples & select indices of clean minorities to oversample
  if (!is.null(r)) {
    n.mino.new <- (r * n - length(ind.mino)) / (1 - r)
    id.mino.new <- rep(ind.mino.clean, times = round(n.mino.new / n.mino.clean))
  } else if (!is.null(dup_size)) {
    n.mino.new <- dup_size * length(ind.mino)
    id.mino.new <- rep(ind.mino.clean, times = floor(n.mino.new / n.mino.clean))
    id.mino.new <- c(id.mino.new,  sample(ind.mino.clean, size = n.mino.new - length(id.mino.new)))
  }
  n.mino.new <- length(id.mino.new)

  # create X for synthetic minority samples
  Xnew <- data.frame(X[id.mino.new, ])
  if (d.num > 0) {
    cons.kernel <- (4 / ((d.num + 2) * n.mino.clean))^(1 / (d.num + 4))
    H <- const * cons.kernel
    Xnew.num <- mvrnorm(n.mino.new, mu = rep(0, d.num), Sigma = covest$cov) * H
    Xnew[, id.num] <- data.matrix(Xnew.num + X[id.mino.new, id.num])
  }

  # create y for synthetic minority samples
  if (classy %in% c("character", "integer", "numeric")) {
    ynew <- as.vector(rep(minoY, NROW(Xnew)), mode = classy)
  } else if (classy == "factor") {
    ynew <- factor(c(rep(minoY, NROW(Xnew))), levels = c(majoY, minoY))
  }

  data.out <- data.frame(ynew, Xnew)
  rownames(data.out) <- 1:NROW(Xnew)

  # re-position columns: put data frame names in the right order
  colnames(data.out) <- colnames(data)[colnames(data) %in% colNames]

  # insert y
  indY <- colnames(data.out) == colNames[1]
  data.out[, indY] <- ynew

  # See whether the order of the variables in formula is the same as in data.
  # If not, then swap columns according to the order in data.
  swap.col <- order(pmatch(colNames[-1], colnames(data.out)[!indY]))
  data.out[, !indY] <- Xnew[, (1:d)[swap.col]]

  out <- list(call            = call,
              ind.mino.out    = ind.mino.out,
              tab.id.mino.new = table(id.mino.new),
              mahdists        = covest$mah,
              data            = data.out)
  class(out) <- "robROSE"
  return(out)
}
