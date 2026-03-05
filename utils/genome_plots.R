#!/usr/bin/env Rscript

### Package management
LoadPackages <- function(req_pkgs, repos = "https://cloud.r-project.org/") {
    for (pkg in req_pkgs) {
        if (!require(pkg, character.only = TRUE)) {
            install.packages(pkg, repos = repos)
            library(pkg, character.only = TRUE)
        }
    }
}

LoadEEMSplotLib <- function() {
    require(devtools)

    if (!require("reemsplots2", character.only = TRUE)) {
        install_github("dipetkov/reemsplots2")
        library("reemsplots2", character.only = TRUE)
    }
}

LoadAsyntLib <- function(commit = NULL) {
    require(stringr)

    tmpdir <- tempdir() # nolint: object_usage_linter.
    system(str_interp("git clone https://github.com/simonhmartin/asynt.git ${tmpdir}"))
    system(str_interp("cd ${tmpdir}; git checkout 1889d10f6598d8b03473531bc396fb08424b9dbe"))
    source(str_interp("${tmpdir}/asynt.R"))
}

### Paths
# plink_path <- "/n/holylfs05/LABS/hopkins_lab/Lab/software/plink_linux/plink"

### Species naming
SpeciesNames <- function(kind = "short", outgroup = FALSE) {
    require(data.table)

    name_tab <- data.table(
        SHORT = c("cusp", "drum", "roem", "pilo"),
        LONG  = c("cuspidata", "drummondii", "roemeriana", "pilosa"),
        ONT   = c("pcusp", "pdrum", "proe", "ppilo"),
        G_SPECIES = c("P. cuspidata", "P. drummondii", "P. roemeriana", "P. pilosa")
    )

    kind <- tolower(kind)
    if (kind %in% c("short", "ddrad")) {
        cn <- "SHORT"
    } else if (kind %in% c("long")) {
        cn <- "LONG"
    } else if (kind %in% c("ont")) {
        cn <- "ONT"
    } else if (kind %in% c("g_species")) {
        cn <- "G_SPECIES"
    } else {
        stop(paste0("kind '", kind, "' not recognized"))
    }

    if (outgroup) {
        out_sp <- name_tab[[cn]]
    } else {
        out_sp <- name_tab[-4, ][[cn]]
    }

    return(out_sp)
}

ConvertSpeciesNames <- function(sp_nm, kind1 = "short", kind2 = "g_species") {
    require(data.table)
    K1 <- K2 <- NULL

    map_tab <- data.table(K1 = SpeciesNames(kind1, outgroup = TRUE),
                          K2 = SpeciesNames(kind2, outgroup = TRUE))
    setkey(map_tab, K1)
    return(map_tab[sp_nm, K2])
}

### Sample naming
RenameWGS <- function(x) {
    #' Rename a long read WGS sample
    require(stringr)

    x <- str_replace(x, "proe_Hopkins_", "R")
    x <- str_replace(x, "pdrum_Hopkins_", "D")
    x <- str_replace(x, "pcusp_Hopkins_", "C")
    x <- str_replace(x, "ppil_Hopkins_", "P")
    x <- str_replace(x, "ppilo_Hopkins_", "P")
    x <- strsplit(x, "_")[[1]][1]
    x <- strsplit(x, "-")[[1]][1]
    x
}

### Data loading
LoadGFF <- function(fn) {
    require(data.table)
    require(stringr)

    if (str_ends(fn, "\\.gz")) {
        use_cmd <- str_interp("gunzip -c ${fn} | grep -v -e '^#' -e '^$'")
    } else {
        use_cmd <- str_interp("grep -v -e '^#' -e '^$' ${fn}")
    }
    out_d <- fread(cmd = use_cmd, header = FALSE, sep = "\t",
                   col.names = c("SEQID", "SRC", "TYPE", "START", "END", "SCORE", "STRAND", "PHASE", "ATTR"))

    out_d
}

LoadPCA <- function(in_prefix) {
    require(data.table)
    require(stringr)

    evec_d <- fread(str_interp("${in_prefix}.eigenvec"), header = FALSE, sep = " ")
    names(evec_d)[1:2] <- c("FID", "IID")
    names(evec_d)[-(1:2)] <- str_c("PC", 1:(ncol(evec_d) - 2))

    eval_d <- fread(str_interp("${in_prefix}.eigenval"), header = FALSE)$V1
    names(eval_d) <- names(evec_d)[-(1:2)]
    prop_var <- eval_d / unname(sum(eval_d))

    list("evec" = evec_d, "eval" = eval_d, "prop_var" = prop_var)
}

LoadTexasBorder <- function() {
    require(rnaturalearth)
    require(sf)

    # Load US states data
    us_states <- ne_states(country = "united states of america", returnclass = "sf")

    # Filter for Texas
    texas <- us_states[us_states[["name"]] == "Texas", ]

    # Extract the coordinates
    texas_coords <- as.data.table(st_coordinates(texas))

    setnames(texas_coords, c("X", "Y"), c("LONG", "LAT"))
    for (cn in str_c("L", 1:3)) {
        texas_coords[, eval(cn) := as.factor(get(cn))]
    }

    texas_coords
}

### BED manipulation
CountUniqueBases <- function(chrom, start_base, end_base, conda_env = "popgen") {
    require(data.table)
    require(stringr)

    tmpfn <- tempfile(fileext = ".bed")
    fwrite(data.table(chrom, start_base, end_base), file = tmpfn, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
    use_cmd <- str_c("conda_setup=$('conda' 'shell.bash' 'hook' 2> /dev/null); eval \"$conda_setup\"; conda activate ", conda_env)
    use_cmd <- str_interp("${use_cmd}; sort -k 1,1 -k2,2n ${tmpfn} | bedtools merge -i - | awk '{sum += $3 - $2} END {print sum}'")
    return(as.numeric(system(use_cmd, intern = TRUE)))
}

### Tree and multiple alignment plotting
ChromEndpoints <- function(chromsizes, gap_size) {
    #' Get the endpoints of the chromosomes
    SIDE <- NULL # for linter
    chromsizes <- as.numeric(chromsizes)
    gap_size <- as.numeric(gap_size)
    x <- cbind(c(0, cumsum(chromsizes)[-length(chromsizes)]), cumsum(chromsizes))
    x <- x + cbind(seq(0, gap_size * (length(chromsizes) - 1), by = gap_size),
                   seq(0, gap_size * (length(chromsizes) - 1), by = gap_size))
    out_list <- data.table(X = as.vector(t(x)),
                           CHROM = as.vector(sapply(str_c("chr", seq_along(chromsizes)), function(x) c(x, x))))
    out_list[, SIDE := rep(c("start", "end"), length(chromsizes))]
    out_list
}

PlotChromosomes <- function(p, ref, width = 0.1, gap_size = 1,
                            label_chrom = FALSE, label_color = "white",
                            label_size = 3, ypos = NULL, ...) {
    #' Plot the chromosome bodies
    require(data.table)
    require(ggplot2)
    SIZE <- SPECIES <- CHROM <- GROUP <- X <- Y <- . <- NULL # for linter
    plt_d <- ref[, ChromEndpoints(SIZE, gap_size), by = SPECIES]
    plt_d <- plt_d[, .(X = c(X, rev(X))), by = .(SPECIES, CHROM)]

    if (is.null(ypos)) {
        plt_d[, Y := as.numeric(as.integer(SPECIES)) + c(-1, -1, 1, 1) * width / 2]
    } else {
        plt_d[, Y := ypos[SPECIES] + c(-1, -1, 1, 1) * width / 2]
    }
    plt_d[, GROUP := str_c(SPECIES, CHROM, sep = "_")]

    p <- p + geom_polygon(aes(x = X, y = Y, group = GROUP, fill = SPECIES), data = plt_d, ...)

    return(p)
}

LabelChromosomes <- function(p, ref, width = 0.1, gap_size = 1, ypos = NULL, ...) {
    #' Label chromosomes with names
    require(data.table)
    X <- Y <- SIZE <- SPECIES <- CHROM <- . <- NULL
    plt_d <- ref[, ChromEndpoints(SIZE, gap_size), by = SPECIES]
    plt_d <- plt_d[, .(X = mean(X)), by = .(SPECIES, CHROM)]

    if (is.null(ypos)) {
        plt_d[, Y := as.numeric(as.integer(SPECIES))]
    } else {
        plt_d[, Y := ypos[SPECIES]]
    }

    p <- p + geom_text(aes(x = X, y = Y, label = CHROM), data = plt_d, ...)
}

BuildConnectorPolygons <- function(reg_dat, end_pos, width, connect_margin, ypos) {
    #' Build curved "riparian" style connector polygons between synteny regions
    require(data.table)

    Rsp <- Rchrom <- Qsp <- Qchrom <- SPECIES <- sp_pair <- NULL # for linter
    Rend <- Rstart  <- Qend <- Qstart <- orient <- NULL # for linter
    idx_a <- reg_dat[, c(Rsp, Rchrom, "start")]
    idx_b <- reg_dat[, c(Qsp, Qchrom, "start")]

    x_offset_a <- end_pos[as.list(idx_a)]$X
    x_offset_b <- end_pos[as.list(idx_b)]$X

    if (is.null(ypos)) {
        y_a <- end_pos[, which(levels(SPECIES) == reg_dat[["Rsp"]])]
        y_b <- end_pos[, which(levels(SPECIES) == reg_dat[["Qsp"]])]
    } else {
        y_a <- end_pos[, ypos[reg_dat[["Rsp"]]]]
        y_b <- end_pos[, ypos[reg_dat[["Qsp"]]]]
    }

    y_offset <- (width / 2) + connect_margin
    if (y_a > y_b) {
        y_offset <- -y_offset
    }

    out_dat <- calc_curvePolygon(
        start1 = reg_dat[["Rstart"]] + x_offset_a,
        end1   = reg_dat[["Rend"]] + x_offset_a,
        start2 = reg_dat[["Qstart"]] + x_offset_b,
        end2   = reg_dat[["Qend"]] + x_offset_b,
        y1     = y_a + y_offset,
        y2     = y_b - y_offset
    )
    out_dat <- as.data.table(out_dat)
    setnames(out_dat, names(out_dat), toupper(names(out_dat)))

    ori <- ifelse(reg_dat[, sign(Rend - Rstart) != sign(Qend - Qstart)], "-", "+")
    out_dat[, sp_pair := reg_dat[["sp_pair"]]]
    out_dat[, orient := ori]

    out_dat
}


PlotConnections <- function(p, dat, ref, width = 0.1, gap_size = 1, connect_margin = 0.01, ypos = NULL, ...) {
    #' Plot riparian style connections
    require(data.table)

    SIZE <- SPECIES <- CHROM <- SIDE <- INDEX <- X <- Y <- orient <- NULL
    end_pos <- ref[, ChromEndpoints(SIZE, gap_size), by = SPECIES]
    setkey(end_pos, SPECIES, CHROM, SIDE)
    plt_d <- dat[, BuildConnectorPolygons(.SD, end_pos, width, connect_margin, ypos), by = INDEX]

    p <- p + geom_polygon(aes(x = X, y = Y, group = INDEX, fill = orient), data = plt_d, ...)
    return(p)
}


### Computation
# QuantileOverReplicates <- function(x, cntr=0.5, ci=c(0.025, 0.975)) {
#     q <- as.list(quantile(x, probs = c(ci[1], cntr, ci[2])))
#     names(q) <- c("LWR", "MED", "UPR")
#     return(q)
# }

QuantileOverReplicates <- function(..., prefixes = "val", sep = "_",
                                   probs = c(0.025, 0.5, 0.975), na.rm = TRUE) { # nolint: object_name_linter.
    in_arg <- list(...)
    if (length(in_arg) != length(prefixes)) {
        stop("Not enough prefixes for provided variables")
    }

    out_list <- list()
    for (i in seq_along(in_arg)) {
        y <- as.list(quantile(in_arg[[i]], probs = probs, na.rm = na.rm))
        names(y) <- str_c(prefixes[i], sep, c("lwr", "med", "upr"))
        out_list <- c(out_list, y)
    }

    return(out_list)
}

BootstrapWtAvg <- function(numer, denom, probs = c(0.025, 0.975), nboot = 1000) {
    # Assumes numer and denom are e.g. measurements across windows.
    # Reports overall average as well as upper and lower percentiles
    # defined in probs
    stopifnot(length(probs) == 2)

    mask <- !is.na(numer) & !is.na(denom)
    numer <- numer[mask]
    denom <- denom[mask]

    n <- length(numer)
    mat <- replicate(sample.int(n, size = n, replace = TRUE), n = nboot)
    boot_avg <- apply(mat, 2, function(x) sum(numer[x]) / sum(denom[x]))
    avg <- sum(numer) / sum(denom)

    q <- unname(quantile(boot_avg, probs = probs))

    list("MEAN" = avg, "LWR" = q[1], "UPR" = q[2])
}

RunOLS <- function(dat, xvar, yvar, npoints = 100, level = 0.95) {
    require(data.table)
    LWR <- UPR <- . <- NULL # for linter

    m <- lm(as.formula(paste(yvar, "~", xvar)), dat)
    pval <- summary(m)$coef[2, 4]
    rsq <- summary(m)$r.squared

    newdat <- dat[, .(X = seq(min(get(xvar)), max(get(xvar)), length.out = npoints))]
    setnames(newdat, "X", xvar)

    tmp <- predict(m, newdat, interval = "confidence", level = level)
    newdat[, eval(yvar) := tmp[, "fit"]]
    newdat[, LWR := tmp[, "lwr"]]
    newdat[, UPR := tmp[, "upr"]]

    list(model = m, pval = pval, rsq = rsq, fitdat = newdat)
}

RunPoisson <- function(dat, xvar, yvar, npoints = 100, level = 0.95,
                       link = "log", ols_start = FALSE, no_intercept = FALSE) {
    require(data.table)
    LWR <- UPR <- . <-  NULL # for linter

    fmla_str <- paste(yvar, "~", xvar)
    if (no_intercept) {
        fmla_str <- paste(fmla_str, " - 1")
    }
    fmla <- as.formula(fmla_str)

    # Run OLS to get starting values
    if (ols_start) {
        ols <- lm(fmla, dat)
        start_vals <- coef(ols)
        start_vals[1] <- max(start_vals[1], .Machine$double.eps)
    } else {
        start_vals <- NULL
    }

    # Fit Poisson GLM
    m <- glm(
        fmla,
        data   = dat,
        family = poisson(link = link),
        start  = start_vals
    )

    # p-value
    pval <- summary(m)$coef[2, 4]

    # Deviance-based pseudo R^2
    rsq <- 1 - m$deviance / m$null.deviance

    newdat <- dat[, .(X = seq(min(get(xvar)), max(get(xvar)), length.out = npoints))]
    setnames(newdat, "X", xvar)

    # Predictions on the link scale, then transform to response
    pr <- predict(m, newdata = newdat, type = "link", se.fit = TRUE)
    linkinv <- m$family$linkinv

    z <- qnorm((1 + level) / 2)
    eta_fit <- pr$fit
    eta_lwr <- eta_fit - z * pr$se.fit
    eta_upr <- eta_fit + z * pr$se.fit

    fit <- linkinv(eta_fit)
    lwr <- linkinv(eta_lwr)
    upr <- linkinv(eta_upr)

    # Fitted values and CIs on the response scale
    newdat[, eval(yvar) := fit]
    newdat[, LWR := lwr]
    newdat[, UPR := upr]

    return(list(model = m, pval = pval, rsq = rsq, fitdat = newdat))
}

RunNegBinom <- function(dat, xvar, yvar, npoints = 100, level = 0.95,
                        link = "log", ols_start = FALSE, no_intercept = FALSE) {
    require(MASS)
    require(data.table)
    LWR <- UPR <- . <- NULL # for linter

    fmla_str <- paste(yvar, "~", xvar)
    if (no_intercept) {
        fmla_str <- paste(fmla_str, " - 1")
    }
    fmla <- as.formula(fmla_str)

    # Run OLS to get starting values
    if (ols_start) {
        ols <- lm(fmla, dat)
        start_vals <- coef(ols)
        start_vals[1] <- max(start_vals[1], .Machine$double.eps)
    } else {
        start_vals <- NULL
    }

    # Fit Negative Binomial GLM
    if (link == "log") {
        m <- MASS::glm.nb(
            formula = fmla,
            data    = dat,
            link    = "log",
            start   = start_vals
        )
    } else if (link == "identity") {
        m <- MASS::glm.nb(
            formula = fmla,
            data    = dat,
            link    = "identity",
            start   = start_vals
        )
    } else if (link == "sqrt") {
        m <- MASS::glm.nb(
            formula = fmla,
            data    = dat,
            link    = "sqrt",
            start   = start_vals
        )
    } else {
        stop("link not recognized")
    }

    pval <- summary(m)$coef[2, 4]

    # Deviance-based pseudo R^2
    rsq <- 1 - m$deviance / m$null.deviance

    newdat <- dat[, .(X = seq(min(get(xvar)), max(get(xvar)), length.out = npoints))]
    data.table::setnames(newdat, "X", xvar)

    # Predictions on the link scale, then transform to response
    pr <- predict(m, newdata = newdat, type = "link", se.fit = TRUE)
    linkinv <- m$family$linkinv

    z <- qnorm((1 + level) / 2)
    eta_fit <- pr$fit
    eta_lwr <- eta_fit - z * pr$se.fit
    eta_upr <- eta_fit + z * pr$se.fit

    fit <- linkinv(eta_fit)
    lwr <- linkinv(eta_lwr)
    upr <- linkinv(eta_upr)

    # Fitted values and CIs on the response scale
    newdat[, (yvar) := fit]
    newdat[, LWR := lwr]
    newdat[, UPR := upr]

    return(list(model = m, pval = pval, rsq = rsq, fitdat = newdat))
}


### Plot stylings

Load4SpPal <- function(pal = "default") {
    if (pal %in% c("default", "bluedrum")) {
        out_colors <- sort(c(roem = "#FFCC00", pilo = "#330066", drum = "#7570B3", cusp = "#339966"))
    } else if (pal %in% c("alternate", "pinkdrum")) {
        out_colors <- sort(c(roem = "#FFCC00", pilo = "#330066", drum = "#993366", cusp = "#339966"))
    } else if (pal %in% c("2drum")) {
        out_colors <- sort(c(roem = "#FFCC00", pilo = "#330066", drumBlue = "#7570B3", drumPink = "#993366", cusp = "#339966"))
    } else {
        stop(str_interp("Palette name ${pal} not recognized"))
    }
    return(out_colors)
}

SetBaseTheme <- function(base_family, fontpath = "fonts/", device = "pdf") {
    require(extrafont)
    require(ggplot2)
    require(ggpubr)

    font_import(paths = fontpath, prompt = FALSE)
    loadfonts(device = device)

    ggplot2::theme_set(theme_pubr(base_family = base_family))
}

GetPlotLimits <- function(p) {
    x <- layer_scales(p)$x$range$range
    y <- layer_scales(p)$y$range$range
    rbind(x, y)
}

### Plotting helper functions
AdmixAnc2Species <- function(plt_d) {
    # Associate STRUCTURE or ADMIXTURE populations to species
    require(data.table)
    V1 <- . <- NULL # for linter

    value <- species <- pop_infer <- NULL # for linter

    pop_infer_totals <- plt_d[, sum(value), by = .(species, pop_infer)]
    admix_pal <- c()
    for (sp in c("roem", "cusp")) {
        pop <- pop_infer_totals[species == sp][which.max(V1), pop_infer]
        admix_pal[pop] <- unname(Load4SpPal("default")[sp])
    }

    remain_pops <- pop_infer_totals[, setdiff(unique(pop_infer), names(admix_pal))]

    list("admix_pal" = admix_pal, "remain_pops" = remain_pops)
}

FloatAsSci <- function(x, dec_lwr = -2, dec_upr = 2, dec_out = 1) {
    # Returns a plotmath string for a value in scientific notation if
    # less than 'dec_lwr'  or greater than 'dec_upr' decimal places.
    # Use 'dec_out' places when outputting in scientific notation
    require(stringr)

    y <- log10(x)
    if (y >= dec_lwr && y <= dec_upr) {
        # Plain text
        if (x >= 1) {
            patt <- str_c("$[", dec_out + 1, ".f]{x}", sep = "")
            outstr <- str_interp(patt)
        } else {
            patt <- str_c("$[0.", dec_out + 1, "f]{x}", sep = "")
            outstr <- str_interp(patt)
        }
    } else {
        # Scientific!
        pwr <- floor(y)
        m <- x * 10^-pwr # nolint: object_usage_linter.
        patt <- str_c("$[0.", dec_out, "f]{m}", sep = "")
        outstr <- c(str_interp(patt), str_interp(" %*% 10^$[d]{pwr}"))
        outstr <- str_c(outstr, collapse = "")
    }
    return(outstr)
}

FormatR2PvalLab <- function(x, stacked = FALSE, pseudo = FALSE) {
    # x should be an [1] rsq value and a [2] pvalue
    require(stringr)

    if (pseudo) {
        r_char <- '"Pseudo-" * italic(R)'
    } else {
        r_char <- "italic(R)"
    }

    if (x[1] == 0) {
        r2_str <- str_c(r_char, "^2==0")
    } else {
        r2_str <- str_c(r_char, "^2==", FloatAsSci(x[1]))
    }

    if (x[2] == 0) {
        pval_str <- str_c("italic(p) < ", FloatAsSci(.Machine$double.eps))
    } else {
        str_c("italic(p)==", FloatAsSci(x[2]))
    }

    if (stacked) {
        out_str <- str_c("atop(", r2_str, ",", pval_str, ")")
    } else {
        out_str <- str_c(r2_str, "~~", pval_str)
    }
    return(out_str)
}

hsv2rgb <- function(x, alpha = FALSE) {
    if (any(is.na(x))) {
        return(rep(NA, 3))
    } else {
        return(grDevices::col2rgb(grDevices::hsv(x[1], x[2], x[3]), alpha = alpha))
    }
}

BrightMatch <- function(x, y) {
    # Match brightness values
    x_hsv <- rgb2hsv(col2rgb(x))
    y_hsv <- rgb2hsv(col2rgb(y))
    x_hsv["v", ] <- y_hsv["v", ]

    tmp <- hsv2rgb(x_hsv) / 255

    rgb(tmp[1], tmp[2], tmp[3])
}

# Metadata2PopulationMap <- function(met_d, min_dist=5) {
#     # Take metadata, identify population centers, and plot

#     pop_d <- met_d[, .(lat=lat[1], long=long[1]), by = Population]
#     pop_dist <- as.data.table(t(combn(pop_d$Population, 2))) # Pairwise distances
#     names(pop_dist) <- c("pop1", "pop2")

#     pop_dist <- merge(pop_dist, pop_d, by.x="pop1", by.y="Population")
#     setnames(pop_dist, c("lat", "long"), c("lat1", "long1"))

#     pop_dist <- merge(pop_dist, pop_d, by.x="pop2", by.y="Population")
#     setnames(pop_dist, c("lat", "long"), c("lat2", "long2"))
#     pop_dist[, DIST := haversine(c(lat1, long1), c(lat2, long2)), by = .(pop1, pop2)]
#     pop_dist[DIST < 1]
# }
