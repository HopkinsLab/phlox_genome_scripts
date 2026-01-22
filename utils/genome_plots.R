#!/usr/bin/env Rscript

### Package management
LoadPackages <- function(req_pkgs, repos="https://cloud.r-project.org/"){
    for(pkg in req_pkgs){
        if(!require(pkg, character.only=TRUE)){install.packages(pkg, repos=repos); library(pkg, character.only=TRUE)}
    }
}

LoadEEMSplotLib <- function(){
    LoadPackages("devtools")
    if(!require("reemsplots2", character.only=TRUE)){
        install_github("dipetkov/reemsplots2")
        library("reemsplots2", character.only=TRUE)
    }
}

### Paths
# plink_path <- "/n/holylfs05/LABS/hopkins_lab/Lab/software/plink_linux/plink"

### Species naming
SpeciesNames <- function(kind="short", outgroup=FALSE){
    LoadPackages("data.table")
    name_tab <- data.table(
        SHORT=c("cusp", "drum", "roem", "pilo"),
        LONG =c("cuspidata", "drummondii", "roemeriana", "pilosa"),
        ONT  =c("pcusp", "pdrum", "proe", "ppilo"),
        G_SPECIES=c("P. cuspidata", "P. drummondii", "P. roemeriana", "P. pilosa")
    )

    kind <- tolower(kind)
    if(kind %in% c("short", "ddrad")){
        cn <- "SHORT"
    } else if(kind %in% c("long")){
        cn <- "LONG"
    } else if(kind %in% c("ont")){
        cn <- "ONT"
    } else if(kind %in% c("g_species")){
        cn <- "G_SPECIES"
    } else {
        stop(paste0("kind '", kind, "' not recognized"))
    }

    if(outgroup){
        out_sp <- name_tab[[cn]]
    } else {
        out_sp <- name_tab[-4,][[cn]]
    }
    
    return(out_sp)
}

ConvertSpeciesNames <- function(sp_nm, kind1="short", kind2="g_species"){
    LoadPackages("data.table")
    map_tab <- data.table(  K1=SpeciesNames(kind1, outgroup=TRUE),
                            K2=SpeciesNames(kind2, outgroup=TRUE))
    setkey(map_tab, K1)
    return(map_tab[sp_nm, K2])
}

### Data loading
LoadGFF <- function(fn){
    LoadPackages(c("data.table", "stringr"))
    if(str_ends(fn, "\\.gz")){
        use_cmd <- str_interp("gunzip -c ${fn} | grep -v -e '^#' -e '^$'")
    } else {
        use_cmd <- str_interp("grep -v -e '^#' -e '^$' ${fn}")
    }
    out_d <- fread(cmd=use_cmd, header=FALSE, sep="\t", 
                col.names=c("SEQID", "SRC", "TYPE", "START", "END", "SCORE", "STRAND", "PHASE", "ATTR"))
    return(out_d)
}

LoadPCA <- function(in_prefix){
    evec_d <- fread(str_interp("${in_prefix}.eigenvec"), header=FALSE, sep=" ")
    names(evec_d)[1:2] <- c("FID", "IID")
    names(evec_d)[-(1:2)] <- str_c("PC", 1:(ncol(evec_d) - 2))

    eval_d <- fread(str_interp("${in_prefix}.eigenval"), header=FALSE)$V1
    names(eval_d) <- names(evec_d)[-(1:2)]
    prop_var <- eval_d / unname(sum(eval_d))

    return(list('evec'=evec_d, 'eval'=eval_d, 'prop_var'=prop_var))
}

LoadTexasBorder <- function(){
    LoadPackages(c("rnaturalearth", "sf"))
    # Load US states data
    us_states <- ne_states(country = "united states of america", returnclass = "sf")

    # Filter for Texas
    texas <- us_states[us_states[["name"]] == "Texas",]
        
    # Extract the coordinates
    texas_coords <- as.data.table(st_coordinates(texas))

    setnames(texas_coords, c("X", "Y"), c("LONG", "LAT"))
    for(cn in str_c("L", 1:3)){
        texas_coords[, eval(cn) := as.factor(get(cn))]
    }
    return(texas_coords)
}

### BED manipulation
CountUniqueBases <- function(chrom, start_base, end_base, conda_env="popgen"){
    LoadPackages(c("data.table", "stringr"))
    tmpfn <- tempfile(fileext=".bed")
    fwrite(data.table(chrom, start_base, end_base), file=tmpfn, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
    use_cmd <- str_c("conda_setup=$('conda' 'shell.bash' 'hook' 2> /dev/null); eval \"$conda_setup\"; conda activate ", conda_env)
    use_cmd <- str_interp("${use_cmd}; sort -k 1,1 -k2,2n ${tmpfn} | bedtools merge -i - | awk '{sum += $3 - $2} END {print sum}'")
    return(as.numeric(system(use_cmd, intern=TRUE)))
}

### Computation
# QuantileOverReplicates <- function(x, cntr=0.5, ci=c(0.025, 0.975)){
#     q <- as.list(quantile(x, probs=c(ci[1], cntr, ci[2])))
#     names(q) <- c("LWR", "MED", "UPR")
#     return(q)
# }

QuantileOverReplicates <- function(..., prefixes="val", sep="_", probs=c(0.025, 0.5, 0.975), na.rm=TRUE){
    in_arg <- list(...)
    if(length(in_arg) != length(prefixes)){
        stop("Not enough prefixes for provided variables")
    }

    out_list <- list()
    for(i in 1:length(in_arg)){
        y <- as.list(quantile(in_arg[[i]], probs=probs, na.rm=na.rm))
        names(y) <- str_c(prefixes[i], sep, c("lwr", "med", "upr"))
        out_list <- c(out_list, y)
    }

    return(out_list)
}

BootstrapWtAvg <- function(numer, denom, probs=c(0.025, 0.975), nboot=1000){
    # Assumes numer and denom are e.g. measurements across windows.
    # Reports overall average as well as upper and lower percentiles
    # defined in probs
    stopifnot(length(probs) == 2)

    mask <- !is.na(numer) & !is.na(denom)
    numer <- numer[mask]
    denom <- denom[mask]

    n <- length(numer)
    mat <- replicate(sample.int(n, size=n, replace=TRUE), n=nboot)
    boot_avg <- apply(mat, 2, function(x) sum(numer[x]) / sum(denom[x]))
    avg <- sum(numer) / sum(denom)

    q <- unname(quantile(boot_avg, probs=probs))
    return(list("MEAN"=avg, "LWR"=q[1], "UPR"=q[2]))
}

RunOLS <- function(dat, xvar, yvar, npoints=100, level=0.95){
    m <- lm(as.formula(paste(yvar, "~", xvar)), dat)
    pval <- summary(m)$coef[2,4]
    rsq <- summary(m)$r.squared

    newdat <- dat[, .(X = seq(min(get(xvar)), max(get(xvar)), length.out=npoints))]
    setnames(newdat, "X", xvar)
    
    tmp <- predict(m, newdat, interval="confidence", level=level)
    newdat[, eval(yvar) := tmp[, "fit"]]
    newdat[, LWR := tmp[, "lwr"]]
    newdat[, UPR := tmp[, "upr"]]

    return(list(model=m, pval=pval, rsq=rsq, fitdat=newdat))
}

### Plot stylings

Load4SpPal <- function(pal="default"){
    if(pal %in% c("default", "bluedrum")){
        out_colors <- sort(c(roem="#FFCC00", pilo="#330066", drum="#7570B3", cusp="#339966"))
    } else if(pal %in% c("alternate", "pinkdrum")){
        out_colors <- sort(c(roem="#FFCC00", pilo="#330066", drum="#993366", cusp="#339966"))
    } else if(pal %in% c("2drum")){
        out_colors <- sort(c(roem="#FFCC00", pilo="#330066", drumBlue="#7570B3", drumPink="#993366", cusp="#339966"))
    } else {
        stop(str_interp("Palette name ${pal} not recognized"))
    }
    return(out_colors)
}

SetBaseTheme <- function(base_family, fontpath="fonts/", device="pdf"){
    LoadPackages(c("extrafont", "ggplot2", "ggpubr"))
    font_import(paths=fontpath, prompt=FALSE)
    loadfonts(device=device)

    ggplot2::theme_set(theme_pubr(base_family=base_family))
}

GetPlotLimits <- function(p){
    x <- layer_scales(p)$x$range$range
    y <- layer_scales(p)$y$range$range
    return(rbind(x, y))
}

### Plotting helper functions
AdmixAnc2Species <- function(plt_d){
    # Associate STRUCTURE or ADMIXTURE populations to species
    LoadPackages("data.table")
    pop_infer_totals <- plt_d[, sum(value), by=.(species, pop_infer)]
    admix_pal <- c()
    for(sp in c("roem", "cusp")){
        pop <- pop_infer_totals[species == sp][which.max(V1), pop_infer]
        admix_pal[pop] <- unname(Load4SpPal("default")[sp])
    }
    
    remain_pops <- pop_infer_totals[, setdiff(unique(pop_infer), names(admix_pal))]

    return(list("admix_pal"=admix_pal, "remain_pops"=remain_pops))
}

FloatAsSci <- function(x, dec_lwr=-2, dec_upr=2, dec_out=1){
    # Returns a plotmath string for a value in scientific notation if
    # less than 'dec_lwr'  or greater than 'dec_upr' decimal places.
    # Use 'dec_out' places when outputting in scientific notation
    LoadPackages("stringr")
    y <- log10(x)
    if(y >= dec_lwr && y <= dec_upr){
        # Plain text
        if(x >= 1){
            patt <- str_c("$[", dec_out+1, ".f]{x}", sep="")
            outstr <- str_interp(patt)
        } else {
            patt <- str_c("$[0.", dec_out+1, "f]{x}", sep="")
            outstr <- str_interp(patt)
        }
    } else {
        # Scientific!
        pwr <- floor(y)
        m <- x * 10^-pwr
        patt <- str_c("$[0.", dec_out, "f]{m}", sep="")
        outstr <- c(str_interp(patt), str_interp(" %*% 10^$[d]{pwr}"))
        outstr <- str_c(outstr, collapse="")
    }
    return(outstr)
}

FormatR2PvalLab <- function(x, stacked=FALSE){
    # x should be an [1] rsq value and a [2] pvalue
    LoadPackages("stringr")
    if(stacked){
        out_str <- str_c("atop(italic(R)^2==", FloatAsSci(x[1]), ",italic(p)==", FloatAsSci(x[2]), ")")
    } else {
        out_str <- str_c("italic(R)^2==", FloatAsSci(x[1]), "~~italic(p)==", FloatAsSci(x[2]))
    }
    return(out_str)
}

hsv2rgb <- function(x, alpha = FALSE) {
  if(any(is.na(x))) {
    return(rep(NA,3)) 
  } else {
    return(grDevices::col2rgb(grDevices::hsv(x[1], x[2], 
                             x[3]), alpha = alpha))
  }
}

BrightMatch <- function(x, y){
    # Match brightness values
    x_hsv <- rgb2hsv(col2rgb(x))
    y_hsv <- rgb2hsv(col2rgb(y))
    x_hsv["v",] <- y_hsv["v",]

    tmp <- hsv2rgb(x_hsv) / 255
    return(rgb(tmp[1], tmp[2], tmp[3]))
}

# Metadata2PopulationMap <- function(met_d, min_dist=5){
#     # Take metadata, identify population centers, and plot

#     pop_d <- met_d[, .(lat=lat[1], long=long[1]), by=Population]
#     pop_dist <- as.data.table(t(combn(pop_d$Population, 2))) # Pairwise distances
#     names(pop_dist) <- c("pop1", "pop2")

#     pop_dist <- merge(pop_dist, pop_d, by.x="pop1", by.y="Population")
#     setnames(pop_dist, c("lat", "long"), c("lat1", "long1"))
    
#     pop_dist <- merge(pop_dist, pop_d, by.x="pop2", by.y="Population")
#     setnames(pop_dist, c("lat", "long"), c("lat2", "long2"))
#     pop_dist[, DIST := haversine(c(lat1, long1), c(lat2, long2)), by=.(pop1, pop2)]
#     pop_dist[DIST < 1]
# }

