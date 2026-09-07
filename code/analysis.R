# =====================================================================
#  Dianchi Lake hydroacoustics - pooled-density reproducible analysis
#  Organised in figure order (Fig 2 -> Fig 6); each block also writes the
#  table(s) that belong with that figure.
#
#  Conversions:  TL = 10^((TS - b0)/m) ;  W = a * TL^b
#  Biomass estimator (pooled method used in the original calculation):
#     mean fish density = arithmetic mean density across valid EDSUs;
#     TS composition = counts pooled across all valid EDSUs;
#     within each length group, a count-weighted mean TS is converted to
#     representative length and weight;
#     whole-lake biomass = mean density x pooled size composition x area.
#
#  Run from a terminal (the working directory does not matter):
#        Rscript path/to/repository/code/analysis.R
#  Or, from the repository root in an interactive R session:
#        source("code/analysis.R", encoding="UTF-8")
#        -> tables in ./tables_analysis/ , figures in ./figures_analysis/
# =====================================================================
script_file <- function() {
  cmd <- grep("^--file=", commandArgs(trailingOnly=FALSE), value=TRUE)
  if (length(cmd)) return(sub("^--file=", "", cmd[[1]]))
  frames <- sys.frames()
  ofiles <- vapply(frames, function(x) {
    if (is.null(x$ofile)) NA_character_ else as.character(x$ofile)
  }, character(1))
  ofiles <- ofiles[!is.na(ofiles)]
  if (length(ofiles)) tail(ofiles, 1) else NA_character_
}

this_script <- script_file()
if (is.na(this_script)) {
  project_dir <- normalizePath(getwd(), winslash="/", mustWork=TRUE)
  if (basename(project_dir) == "code") project_dir <- dirname(project_dir)
} else {
  project_dir <- dirname(dirname(normalizePath(this_script, winslash="/", mustWork=TRUE)))
}

need <- c("readxl","ggplot2","patchwork","scales","MCMCglmm","ggrepel")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly=TRUE)]
if (length(miss)) {
  stop(sprintf("Missing required R package(s): %s. See README.md for installation instructions.",
               paste(miss, collapse=", ")), call.=FALSE)
}
invisible(lapply(need, function(p) suppressMessages(library(p, character.only=TRUE))))
base   <- project_dir
DATA   <- file.path(base,"data")
outdir <- file.path(base,"tables_analysis");  dir.create(outdir, showWarnings=FALSE)
figdir <- file.path(base,"figures_analysis"); dir.create(figdir, showWarnings=FALSE)
FIG2_FILE <- file.path(DATA, "fish_catch_data.xlsx")
if (!file.exists(FIG2_FILE)) stop("Cannot find data/fish_catch_data.xlsx.")
set.seed(2025)

A_survey_ha <- 33000     # surveyed lake area = 330 km2
N_SITES     <- 12
BOOT_N      <- 2000; BLOCK <- 20
BIO_LABELS  <- c("4.8-10.2","10.2-21.9","21.9-36.3","36.3-60.3","60.3-129.2")
TS_HEADER_LABELS <- c(-64,-61,-58,-56,-53,-50,-48,-45,-43,-40,-37,-35,-32,-29)
# The header workbook displays rounded TS labels. These Foote-derived class
# centres reproduce the more precise midpoints used in the original spreadsheet.
TS_CLASS_TL_FOOTE <- c(2.58,3.49,4.74,6.42,8.70,11.80,15.99,
                       21.68,29.39,39.84,54.01,73.22,99.26,134.56)
TS_GROUP_ID <- c(1,1,1,1, 2,2,2, 3,3, 4,4, 5,5,5)

## ---- TS-length conversion formulas (Table 2) -----------------------
TSL <- data.frame(name=c("Chen et al. (2019)","Frouzova et al. (2005)","Love (1971)",
                         "Foote (1987)","Ren et al. (2011)"),
                  m=c(35.88,23.97,19.10,20.00,23.50), b0=c(-90.33,-79.93,-63.90,-71.90,-87.20),
                  freq=c(120,120,120,38,200), stringsAsFactors=FALSE)
sel  <- match("Frouzova et al. (2005)", TSL$name)   # equation selected in Section 3.4
TLen <- function(m,b0,ts) 10^((ts-b0)/m)
Wclass <- function(m,b0,a,b) a*TLen(m,b0,ts_mid)^b  # per-class weight (g)

## Species-name lookup (source records use Chinese names; outputs use Latin names).
sp_map <- tryCatch(read.csv(file.path(DATA,"species_names.csv"),
                            encoding="UTF-8", stringsAsFactors=FALSE), error=function(e) NULL)
lat <- function(x){ if(is.null(sp_map)) return(x)
  i <- match(x, sp_map$chinese); ifelse(is.na(i), x, sp_map$latin[i]) }

## ---- acoustic EDSU data: per-cell TS-class counts + density --------
ac  <- suppressMessages(read_excel(file.path(DATA,"hydroacoustic_edsu_data.xlsx")))
nm  <- suppressWarnings(as.numeric(names(ac)))
cls <- which(!is.na(nm) & nm > -70 & nm < -20)
if (length(cls) != 14) stop(sprintf("Expected 14 TS columns, found %d.", length(cls)))
cls <- cls[order(nm[cls])]
ts_header <- nm[cls]
ts_mid <- ts_header
if (isTRUE(all.equal(ts_header, TS_HEADER_LABELS, tolerance=1e-8))) {
  ts_mid <- 20*log10(TS_CLASS_TL_FOOTE)-71.90
  message("Rounded TS headers detected; using calibrated class midpoints from the original calculation.")
}
density_name <- intToUtf8(c(23494, 24230))
if (density_name %in% names(ac)) {
  dcol <- density_name
} else {
  density_candidates <- grep("density", names(ac), ignore.case=TRUE)
  if (!length(density_candidates)) {
    stop("Cannot identify the fish-density column in the hydroacoustic workbook.")
  }
  dcol <- names(ac)[density_candidates[[1]]]
}
den <- as.numeric(ac[[dcol]]); Cnt <- as.matrix(ac[,cls]); Cnt[is.na(Cnt)] <- 0
keep<- rowSums(Cnt)>0 & !is.na(den); Cnt <- Cnt[keep,,drop=FALSE]; den <- den[keep]
nE  <- nrow(Cnt); cls_tot <- colSums(Cnt)
cat(sprintf("Valid EDSUs: %d ; TS classes: %d\n", nE, length(cls)))

## ---- complete catch workbook and measured-individual subset --------
## All sheets are sampling sites. Rows with N > 1 are aggregate catch
## records: they contribute to Fig. 2a-b and IRI, but never to individual
## length-weight analyses. Fig. 2c onward uses exactly the same strict
## subset: N = 1 and positive TL, SL and W.
catch_sheets <- excel_sheets(FIG2_FILE)
if (length(catch_sheets) != N_SITES)
  stop(sprintf("Expected %d catch-site sheets, found %d.", N_SITES, length(catch_sheets)))
catch_raw <- do.call(rbind, lapply(catch_sheets, function(sh){
  d <- read_excel(FIG2_FILE, sheet=sh)
  if (ncol(d) < 5) stop(sprintf("Catch sheet '%s' has fewer than five columns.", sh))
  d <- as.data.frame(d[,1:5], stringsAsFactors=FALSE)
  names(d) <- c("species","TL","SL","W","N")
  d$site <- sh
  d
}))
catch_raw$species <- trimws(as.character(catch_raw$species))
for (z in c("TL","SL","W","N"))
  catch_raw[[z]] <- suppressWarnings(as.numeric(catch_raw[[z]]))

## Complete catch totals for Fig. 2a-b and IRI.
fig2_all <- catch_raw[!is.na(catch_raw$species) & nzchar(catch_raw$species) &
                      is.finite(catch_raw$W) & catch_raw$W > 0 &
                      is.finite(catch_raw$N) & catch_raw$N > 0,]

## Strict complete-case individual data used consistently from Fig. 2c onward.
fig2_ind <- fig2_all[fig2_all$N == 1 &
                     is.finite(fig2_all$TL) & fig2_all$TL > 0 &
                     is.finite(fig2_all$SL) & fig2_all$SL > 0 &
                     is.finite(fig2_all$W)  & fig2_all$W  > 0,]
if (!nrow(fig2_ind)) stop("No complete individual records (N=1, TL/SL/W > 0) were found.")
fig2_lw <- fig2_all[fig2_all$N == 1 &
                    is.finite(fig2_all$TL) & fig2_all$TL > 0 &
                    is.finite(fig2_all$W)  & fig2_all$W  > 0,]
filter_audit <- data.frame(
  item=c("site_sheets","raw_rows","valid_rows_for_Fig2a_b",
         "total_fish_for_Fig2a_b","complete_TL_W_individuals",
         "complete_TL_SL_W_individuals","excluded_for_missing_SL"),
  value=c(length(catch_sheets),nrow(catch_raw),nrow(fig2_all),
          sum(fig2_all$N),nrow(fig2_lw),nrow(fig2_ind),
          nrow(fig2_lw)-nrow(fig2_ind)),
  stringsAsFactors=FALSE)
write.csv(filter_audit, file.path(outdir,"TableS_data_filter_audit.csv"),
          row.names=FALSE, fileEncoding="UTF-8")

cc <- fig2_ind
cc$TLn <- cc$TL
cc$Wn  <- cc$W
gl_len <- cc$TLn
CATCH_MAX <- max(gl_len, na.rm=TRUE)

## The former fitted-data workbook was only a two-column copy of the old catch.
## Refit the pooled community relationship directly from the new complete
## measured-individual dataset so every downstream use has one data source.
fit <- data.frame(TL=cc$TLn, W=cc$Wn)
lp  <- lm(log10(W)~log10(TL), fit)
community_a <- unname(10^coef(lp)[1])
community_b <- unname(coef(lp)[2])
community_r2 <- summary(lp)$r.squared

## ---- length-weight conversion formulas (Table 1) -------------------
LWR <- data.frame(name=c("Community-fitted","Ye et al. (2007)","Wanner and Klumb (2009)"),
                  a=c(community_a,0.0052,0.013),
                  b=c(community_b,3.162,2.96), stringsAsFactors=FALSE)
write.csv(LWR, file.path(outdir,"Table1_length_weight_equations.csv"),
          row.names=FALSE, fileEncoding="UTF-8")

## ---- pooled biomass helper (mean density x pooled TS composition) --
pooled_group_bio <- function(m,b0,a,b, d=den, C=Cnt, return_groups=FALSE){
  counts <- colSums(C, na.rm=TRUE)
  total_counts <- sum(counts)
  mean_density <- mean(d, na.rm=TRUE)
  if (!is.finite(mean_density) || total_counts <= 0)
    stop("Cannot calculate pooled biomass: density or TS counts are invalid.")

  # Fixed 4+3+2+2+3 TS-class grouping, matching the original spreadsheet.
  class_group <- TS_GROUP_ID
  ng <- length(BIO_LABELS)
  group_counts <- vapply(seq_len(ng), function(g)
    sum(counts[class_group==g]), numeric(1))
  group_mean_ts <- vapply(seq_len(ng), function(g){
    z <- class_group==g
    if (!any(z) || sum(counts[z]) <= 0) return(NA_real_)
    weighted.mean(ts_mid[z], counts[z])
  }, numeric(1))

  # This follows the spreadsheet: group mean TS -> representative TL -> W.
  group_weight_g <- a * TLen(m,b0,group_mean_ts)^b
  group_density <- mean_density * group_counts / total_counts
  group_tonnes <- group_density * group_weight_g * A_survey_ha / 1e6
  group_tonnes[!is.finite(group_tonnes)] <- 0
  if (return_groups) group_tonnes else sum(group_tonnes)
}

## palettes / theme (match the manuscript figures) -------------------
tsl_short <- c("Chen 2019","Frouzova 2005","Love 1971","Foote 1987","Ren 2011")
lwr_short <- c("Community-fitted","Ye 2007","Wanner & Klumb 2009")
col_num <- "#3C5488"; col_bio <- "#E64B35"
pal_lwr <- setNames(c("#3C5488","#E64B35","#E1A000"), lwr_short)
pal_tsl <- setNames(c("#E64B35","#3C5488","#00A087","#F39B7F","#8491B4"), tsl_short)
gld_lev <- c("Planktivore","Carnivore","Omnivore","Filter-feeder")
pal_gld <- setNames(c("#92C5DE","#2166AC","#F4A582","#B2182B"), gld_lev)
th <- theme_bw(base_size=10) +
      theme(panel.grid.minor=element_blank(), plot.tag=element_text(face="bold"),
            legend.title=element_text(size=8), legend.text=element_text(size=7),
            legend.key.size=unit(0.4,"cm"))
sav <- function(p, f, w, h) {
  
  stem <- tools::file_path_sans_ext(f)
  
  # PNG files are convenient for routine viewing.
  ggsave(
    filename = file.path(figdir, paste0(stem, ".png")),
    plot = p,
    width = w,
    height = h,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  
  # High-resolution TIFF files are suitable for manuscript submission.
  ggsave(
    filename = file.path(figdir, paste0(stem, ".tiff")),
    plot = p,
    device = "tiff",
    width = w,
    height = h,
    units = "in",
    dpi = 600,
    compression = "lzw",
    bg = "white"
  )
}
                                  
# =====================================================================
#  FIG 2 - community structure from the updated catch workbook
#  Panels a-b: all valid catch totals (including aggregate-only rows).
#  Panels c-d: individually measured fish only (N=1, TL/SL/W complete).
#  The same complete individual subset `cc` is used from Fig 3 onward.
# =====================================================================
t3 <- read.csv(file.path(DATA,"table3_dominance.csv"), fileEncoding="UTF-8", stringsAsFactors=FALSE)

map_fig2 <- function(d){
  d$latin <- lat(d$species)
  d$code <- t3$code[match(d$latin, t3$species)]
  d$feeding_type <- t3$feeding_type[match(d$code, t3$code)]
  missing_species <- sort(unique(d$species[is.na(d$code) | is.na(d$feeding_type)]))
  if (length(missing_species))
    stop("Missing species/guild mapping for Fig. 2: ",
         paste(missing_species, collapse=", "))
  d
}
fig2_all <- map_fig2(fig2_all)
fig2_ind <- map_fig2(fig2_ind)
fig2_ind$guild <- factor(fig2_ind$feeding_type, levels=gld_lev)
cc <- fig2_ind
cc$TLn <- cc$TL
cc$Wn  <- cc$W
cc$latin <- lat(cc$species)
cc$code <- as.character(cc$code)
cc$guild <- factor(cc$feeding_type, levels=gld_lev)

## Complete-catch dominance for Fig. 2a-b and the updated Table 3.
dom_all <- aggregate(cbind(N,W)~code+latin+feeding_type, fig2_all, sum)
dom_all$N_pct <- dom_all$N / sum(dom_all$N) * 100
dom_all$W_pct <- dom_all$W / sum(dom_all$W) * 100
occ <- aggregate(site~code, unique(fig2_all[,c("code","site")]), length)
names(occ)[2] <- "sites_present"
dom_all <- merge(dom_all, occ, by="code", all.x=TRUE, sort=FALSE)
dom_all$F <- dom_all$sites_present / N_SITES
dom_all$IRI <- (dom_all$N_pct/100 + dom_all$W_pct/100) * dom_all$F * 10000
dom_all$dominant <- ifelse(dom_all$IRI >= 1000, "yes", "no")
dom_all <- dom_all[order(-dom_all$N_pct),]
write.csv(dom_all[,c("latin","code","feeding_type","N","W","N_pct","W_pct",
                    "sites_present","F","IRI","dominant")],
          file.path(outdir,"Table3_updated_catch_dominance.csv"),
          row.names=FALSE, fileEncoding="UTF-8")
cat(sprintf(paste0("Fig. 2 updated catch: %d fish, %.1f g; ",
                   "complete individually measured subset: %d fish; catch maximum TL: %.1f cm.\n"),
            sum(fig2_all$N), sum(fig2_all$W), nrow(fig2_ind), CATCH_MAX))

## Fig 2a : dominance decoupling (complete catch, ordered by numerical %)
lev_num <- dom_all$code[order(dom_all$N_pct)]
d2a <- rbind(data.frame(code=dom_all$code, pct=dom_all$N_pct, metric="Numerical %"),
             data.frame(code=dom_all$code, pct=dom_all$W_pct, metric="Biomass %"))
d2a$code   <- factor(d2a$code, levels=lev_num)
d2a$metric <- factor(d2a$metric, levels=c("Numerical %","Biomass %"))
p2a <- ggplot() +
  geom_segment(data=dom_all,
               aes(y=factor(code,lev_num), yend=factor(code,lev_num),
                   x=N_pct, xend=W_pct),
               colour="grey70", linewidth=0.5) +
  geom_point(data=d2a, aes(y=code, x=pct, colour=metric), size=2.2) +
  scale_colour_manual(values=c("Numerical %"=col_num,"Biomass %"=col_bio)) +
  labs(x="Percentage (%)", y=NULL, colour=NULL, tag="a") +
  th + theme(legend.position="inside", legend.position.inside=c(0.8,0.35))

## Fig 2b : feeding-guild composition (complete catch)
gsum <- aggregate(cbind(N_pct,W_pct)~feeding_type, dom_all, sum)
d2b <- rbind(data.frame(g=gsum$feeding_type, pct=gsum$N_pct, metric="Numerical %"),
             data.frame(g=gsum$feeding_type, pct=gsum$W_pct, metric="Biomass %"))
d2b$g <- factor(d2b$g, levels=rev(gld_lev))
d2b$metric <- factor(d2b$metric, levels=c("Numerical %","Biomass %"))
p2b <- ggplot(d2b, aes(pct, g, fill=metric)) +
  geom_col(position=position_dodge(0.7), width=0.65) +
  scale_fill_manual(values=c("Numerical %"=col_num,"Biomass %"=col_bio)) +
  labs(x="Percentage of total", y=NULL, fill=NULL, tag="b") +
  th + theme(legend.position="inside", legend.position.inside=c(0.78,0.55))

## Fig 2c : dominance vs mean body size in the individually measured subset
dom_ind <- aggregate(cbind(N,W)~code, fig2_ind, sum)
mean_tl_ind <- aggregate(TL~code, fig2_ind, mean)
dom_ind <- merge(dom_ind, mean_tl_ind, by="code", all=FALSE)
dom_ind$N_pct <- dom_ind$N / sum(dom_ind$N) * 100
dom_ind$W_pct <- dom_ind$W / sum(dom_ind$W) * 100
d2c <- rbind(data.frame(code=dom_ind$code, mtl=dom_ind$TL,
                        pct=dom_ind$N_pct, metric="Numerical %"),
             data.frame(code=dom_ind$code, mtl=dom_ind$TL,
                        pct=dom_ind$W_pct, metric="Biomass %"))
d2c <- d2c[is.finite(d2c$mtl) & d2c$pct>0,]
d2c$metric <- factor(d2c$metric, levels=c("Numerical %","Biomass %"))
p2c <- ggplot(d2c, aes(mtl, pct, colour=metric, fill=metric)) +
  geom_smooth(method="lm", formula=y~x, se=TRUE, linewidth=0.8, alpha=0.15) +
  geom_point(size=1.8) +
  geom_text_repel(aes(label=code), size=2.3, show.legend=FALSE,
                  max.overlaps=30, seed=1) +
  scale_colour_manual(values=c("Numerical %"=col_num,"Biomass %"=col_bio)) +
  scale_fill_manual(values=c("Numerical %"=col_num,"Biomass %"=col_bio)) +
  scale_x_log10() + scale_y_log10() +
  labs(x="Mean total length (cm)", y="Relative dominance (%)", tag="c") +
  th + theme(legend.position="none")
## dominance-size regression slopes (reported in Section 3.1)
dsd <- dom_ind[is.finite(dom_ind$TL) & dom_ind$N_pct>0 & dom_ind$W_pct>0,]
snum <- lm(log10(N_pct)~log10(TL), dsd)
sbio <- lm(log10(W_pct)~log10(TL), dsd)
cross_cm <- 10^((coef(snum)[1]-coef(sbio)[1])/(coef(sbio)[2]-coef(snum)[2]))

## Fig 2d : biomass composition by length group and guild, measured fish only
bnd2  <- c(2.2,7.9,15.7,29.4,54.0,200)
glab2 <- c("2.2-7.9","7.9-15.7","15.7-29.4","29.4-54.0","54.0-134.6")
fig2_ind$lg <- cut(fig2_ind$TL, breaks=bnd2, labels=glab2,
                   include.lowest=TRUE, right=FALSE)
d2d <- aggregate(W~lg+guild, fig2_ind[!is.na(fig2_ind$guild),], sum)
d2d$frac <- d2d$W / as.numeric(tapply(d2d$W, d2d$lg, sum)[as.character(d2d$lg)])
d2d$lg <- factor(d2d$lg, levels=rev(glab2))
d2d$guild <- factor(d2d$guild, levels=gld_lev)
p2d <- ggplot(d2d, aes(frac, lg, fill=guild)) +
  geom_col(width=0.72, colour="white", linewidth=0.2) +
  scale_fill_manual(values=pal_gld) + scale_x_continuous(labels=label_percent()) +
  labs(x="Biomass composition", y="Length group (cm)", fill=NULL, tag="d") +
  th + theme(legend.position="bottom")
sav((p2a|p2b)/(p2c|p2d), "Fig2_community_structure.png", 10, 8)

# =====================================================================
#  FIG 3 - length-weight relationships (Tables 4 & 5)
# =====================================================================
cat(sprintf("Pooled LWR: a=%.5f b=%.3f R2=%.3f (n=%d)\n",
            community_a, community_b, community_r2, nrow(fit)))

## Table 5 : species-specific LWR (n >= 20)
res <- do.call(rbind, lapply(split(cc, cc$code), function(d){
  if(nrow(d)<20) return(NULL); m<-lm(log10(Wn)~log10(TLn),d); ci<-confint(m)[2,]
  data.frame(species=lat(d$species[1]), n=nrow(d), a=round(10^coef(m)[1],5), b=round(coef(m)[2],3),
             b_95CI=sprintf("%.3f-%.3f",ci[1],ci[2]), R2=round(summary(m)$r.squared,3),
             type=ifelse(ci[2]<3,"-",ifelse(ci[1]>3,"+","I")))}))
res <- res[order(res$b),]
write.csv(res, file.path(outdir,"Table5_species_LWR.csv"), row.names=FALSE, fileEncoding="UTF-8")

## Heterogeneity-of-slopes ANCOVA for species with n >= 20.
sp_n <- table(cc$code)
cc20 <- cc[cc$code %in% names(sp_n[sp_n >= 20]),]
cc20$code <- droplevels(factor(cc20$code))
anc_lwr <- lm(log10(Wn) ~ log10(TLn) * code, data=cc20)
anc_tab <- anova(anc_lwr)
anc_row <- "log10(TLn):code"
anc_F <- anc_tab[anc_row,"F value"]
anc_df <- anc_tab[anc_row,"Df"]
anc_resid_df <- anc_tab["Residuals","Df"]
anc_P <- anc_tab[anc_row,"Pr(>F)"]
anc_out <- data.frame(
  test="log10(W) ~ log10(TL) x species",
  species_n_min=20,
  interaction_df=anc_df,
  residual_df=anc_resid_df,
  F=anc_F,
  P=anc_P,
  stringsAsFactors=FALSE)
write.csv(anc_out, file.path(outdir,"TableS_species_slope_heterogeneity.csv"),
          row.names=FALSE, fileEncoding="UTF-8")
cat(sprintf("Species-slope ANCOVA: F=%.2f, df=%d and %d, P=%s\n",
            anc_F, anc_df, anc_resid_df,
            ifelse(anc_P<0.001,"<0.001",sprintf("%.4f",anc_P))))

## Table 4 : length and weight statistics
t4 <- do.call(rbind, lapply(split(cc, cc$code), function(d) data.frame(
  species=lat(d$species[1]), n=nrow(d),
  TL_range_cm=sprintf("%.1f-%.1f", min(d$TLn), max(d$TLn)),
  mean_TL_cm=sprintf("%.2f +/- %.2f", mean(d$TLn), sd(d$TLn)),
  weight_range_g=sprintf("%.1f-%.0f", min(d$Wn), max(d$Wn)),
  mean_weight_g=sprintf("%.1f +/- %.1f", mean(d$Wn), sd(d$Wn)), stringsAsFactors=FALSE)))
t4 <- t4[order(-t4$n), ]
write.csv(t4, file.path(outdir,"Table4_TL_weight_stats.csv"), row.names=FALSE, fileEncoding="UTF-8")

## hierarchical Bayesian length-weight exponent (MCMCglmm; Hadfield 2010)
## log10(W) ~ log10(TL) with species-varying intercept and slope, partially
## pooled toward a community mean; weakly informative parameter-expanded priors.
cc_m <- cc[!is.na(cc$code),]
cc_m$logW <- log10(cc_m$Wn); cc_m$logTL <- log10(cc_m$TLn); cc_m$sp <- factor(cc_m$code)
bprior <- list(R=list(V=1, nu=0.002),
               G=list(G1=list(V=diag(2), nu=2, alpha.mu=rep(0,2), alpha.V=diag(2)*1000)))
set.seed(2025)
mm <- MCMCglmm(logW~logTL, random=~us(1+logTL):sp, data=cc_m, prior=bprior,
               nitt=23000, burnin=3000, thin=20, pr=TRUE, verbose=FALSE)
comm_b  <- mean(mm$Sol[,"logTL"]); comm_ci <- as.numeric(HPDinterval(mm$Sol[,"logTL"]))
fixsl <- mm$Sol[,"logTL"]                              # posterior of the population slope
bayes <- do.call(rbind, lapply(levels(cc_m$sp), function(cd){   # per-species posterior slopes
  rc <- paste0("logTL.sp.",cd); if(!(rc %in% colnames(mm$Sol))) return(NULL)
  ps <- fixsl + mm$Sol[,rc]
  data.frame(code=cd, b=mean(ps), lo=HPDinterval(ps)[1], hi=HPDinterval(ps)[2], stringsAsFactors=FALSE)}))
cat(sprintf("Bayesian community exponent %.3f (95%% CrI %.3f-%.3f)\n", comm_b, comm_ci[1], comm_ci[2]))

## length-weight accuracy of each formula against the measured catch (Section 3.2)
acc1 <- function(a,b){ pr<-a*cc$TLn^b
  c(total_bias_pct=(sum(pr)/sum(cc$Wn)-1)*100, median_err_pct=median((pr-cc$Wn)/cc$Wn*100),
    log10_RMSE=sqrt(mean((log10(pr)-log10(cc$Wn))^2))) }
accT <- rbind(Community=acc1(LWR$a[1],LWR$b[1]), Ye_2007=acc1(LWR$a[2],LWR$b[2]), Wanner_2009=acc1(LWR$a[3],LWR$b[3]))
# 10-fold cross-validation of the community fit (out-of-sample, avoids circularity)
set.seed(2025)
folds <- sample(rep(1:10, length.out=nrow(cc))); cvpred <- numeric(nrow(cc))
for(k in 1:10){ tr<-folds!=k; m<-lm(log10(Wn)~log10(TLn),cc[tr,]); cvpred[!tr]<-10^predict(m,cc[!tr,]) }
cv_rmse <- sqrt(mean((log10(cvpred)-log10(cc$Wn))^2)); cv_bias <- (sum(cvpred)/sum(cc$Wn)-1)*100
accTab <- data.frame(formula=c(rownames(accT),"Community (10-fold CV)"),
  total_bias_pct=round(c(accT[,1],cv_bias),1), median_err_pct=round(c(accT[,2],NA),1),
  log10_RMSE=round(c(accT[,3],cv_rmse),3))
write.csv(accTab, file.path(outdir,"TableS_LWR_accuracy.csv"), row.names=FALSE, fileEncoding="UTF-8")

## Fig 3a : pooled fit + 95% band (linear axes)
lg <- data.frame(TL=seq(min(fit$TL), max(fit$TL), length.out=200))
pr <- predict(lp, newdata=lg, interval="confidence")
lg$W <- 10^pr[,"fit"]; lg$lo <- 10^pr[,"lwr"]; lg$hi <- 10^pr[,"upr"]
p3a <- ggplot() +
  geom_point(data=fit, aes(TL,W), colour="grey65", size=0.5, alpha=0.5) +
  geom_ribbon(data=lg, aes(TL, ymin=lo, ymax=hi), fill=col_bio, alpha=0.35) +
  geom_line(data=lg, aes(TL,W), colour=col_bio, linewidth=0.9) +
  annotate("text", x=4, y=max(fit$W)*0.98, hjust=0, vjust=1, size=3.4, parse=TRUE,
           label=sprintf("italic(W)==%.5f*italic(L)^{%.3f}",community_a,community_b)) +
  annotate("text", x=4, y=max(fit$W)*0.80, hjust=0, vjust=1, size=3.4, parse=TRUE,
           label=sprintf("italic(R)^{2}==%.3f",community_r2)) +
  scale_y_continuous(labels=label_comma()) +
  labs(x="Total length (cm)", y="Body weight (g)", tag="a") + th

## Fig 3b : three formulas (linear axes, distinct line types)
lx <- seq(2,135,length.out=300)
d3b <- do.call(rbind, Map(function(a,b,nm) data.frame(TL=lx, W=a*lx^b, f=nm), LWR$a, LWR$b, lwr_short))
d3b$f <- factor(d3b$f, levels=lwr_short)
p3b <- ggplot(d3b, aes(TL,W,colour=f,linetype=f)) + geom_line(linewidth=0.9) +
  scale_colour_manual(values=pal_lwr) +
  scale_linetype_manual(values=setNames(c("solid","dashed","dotdash"),lwr_short)) +
  scale_y_continuous(labels=label_comma()) +
  labs(x="Total length (cm)", y="Predicted weight (g)", colour=NULL, linetype=NULL, tag="b") +
  th + theme(legend.position="inside", legend.position.inside=c(0.32,0.85))

## Fig 3c : species lines (log-log, no legend)
spd <- split(cc, cc$code); spd <- spd[vapply(spd,nrow,0L)>=20]
d3c <- do.call(rbind, lapply(spd, function(d){
  m<-lm(log10(Wn)~log10(TLn),d); xx<-exp(seq(log(min(d$TLn)),log(max(d$TLn)),length.out=50))
  data.frame(TL=xx, W=10^(coef(m)[1]+coef(m)[2]*log10(xx)), code=d$code[1])}))
p3c <- ggplot() +
  geom_point(data=cc, aes(TLn,Wn), colour="grey75", size=0.4, alpha=0.4) +
  geom_line(data=d3c, aes(TL,W,colour=code), linewidth=0.7) +
  scale_x_log10() + scale_y_log10(labels=label_comma()) +
  labs(x="Total length (cm, log)", y="Body weight (g, log)", tag="c") +
  th + theme(legend.position="none")

## Fig 3d : OLS vs hierarchical-Bayesian exponent
ols <- do.call(rbind, lapply(split(cc_m,cc_m$code), function(d){
  if(nrow(d)<3) return(NULL); m<-lm(log10(Wn)~log10(TLn),d); ci<-confint(m)[2,]
  data.frame(code=d$code[1], b=coef(m)[2], lo=ci[1], hi=ci[2])}))
d3d <- merge(ols, bayes, by="code", suffixes=c("_o","_b")); d3d <- d3d[order(d3d$b_o),]; codelev <- d3d$code
d3d_l <- rbind(data.frame(code=d3d$code, b=d3d$b_o, lo=d3d$lo_o, hi=d3d$hi_o, method="OLS"),
               data.frame(code=d3d$code, b=d3d$b_b, lo=d3d$lo_b, hi=d3d$hi_b, method="Bayesian"))
d3d_l$method <- factor(d3d_l$method, levels=c("OLS","Bayesian"))
d3d_l$yy <- match(d3d_l$code, codelev) + ifelse(d3d_l$method=="OLS", 0.16, -0.16)
p3d <- ggplot(d3d_l) +
  geom_vline(xintercept=3, linetype=3, colour="grey40") +
  geom_vline(xintercept=coef(lp)[2], linetype=2, colour=col_num) +
  geom_vline(xintercept=comm_b, linetype=1, colour=col_bio) +
  geom_segment(aes(x=lo, xend=hi, y=yy, yend=yy, colour=method), linewidth=0.45) +
  geom_point(aes(x=b, y=yy, colour=method, shape=method), size=2.2) +
  scale_colour_manual(values=c("OLS"=col_num,"Bayesian"=col_bio), name=NULL) +
  scale_shape_manual(values=c("OLS"=16,"Bayesian"=18), name=NULL) +
  scale_y_continuous(breaks=seq_along(codelev), labels=codelev) +
  coord_cartesian(xlim=c(1.9,4.1)) +
  annotate("text", x=comm_b, y=length(codelev)+0.4,
           label=sprintf("community %.2f",comm_b), colour=col_bio, size=2.5, hjust=1.05) +
  annotate("text", x=coef(lp)[2], y=0.7,
           label=sprintf("pooled %.2f",coef(lp)[2]), colour=col_num, size=2.5, hjust=-0.05) +
  labs(x="Growth exponent b (+/-95% interval)", y=NULL, tag="d") +
  th + theme(legend.position="inside", legend.position.inside=c(0.86,0.16),
             legend.background=element_rect(fill="white", colour="grey80", linewidth=0.2),
             legend.margin=margin(1,3,1,3), axis.text.y=element_text(size=8))
sav((p3a|p3b)/(p3c|p3d), "Fig3_length_weight.png", 10, 8)

# =====================================================================
#  FIG 4 - consistency-based selection of the TS-length equation (Table 2)
# =====================================================================
maxTL <- TLen(TSL$m, TSL$b0, max(ts_mid))
ks_D  <- sapply(seq_len(nrow(TSL)), function(i){
  L <- TLen(TSL$m[i],TSL$b0[i],ts_mid); ac_len <- rep(L, times=round(cls_tot))
  suppressWarnings(ks.test(gl_len, ac_len)$statistic) })
seltab <- data.frame(TS_length_equation=TSL$name, frequency_kHz=TSL$freq,
                     max_inferred_TL_cm=round(maxTL),
                     reproduces_catch_max=ifelse(maxTL>=CATCH_MAX,"yes","no"),
                     KS_distance_D=round(ks_D,3))
write.csv(seltab, file.path(outdir,"Fig4_TSlength_selection.csv"), row.names=FALSE, fileEncoding="UTF-8")
cat("K-S distances:", paste(sprintf("%s=%.2f",TSL$name,ks_D),collapse="; "), "\n")

tlx <- exp(seq(log(2), log(135), length.out=200))
d4a <- do.call(rbind, lapply(seq_len(nrow(TSL)), function(i)
  data.frame(TL=tlx, TS=TSL$m[i]*log10(tlx)+TSL$b0[i], eq=tsl_short[i])))
d4a$eq <- factor(d4a$eq, levels=tsl_short)
p4a <- ggplot(d4a, aes(TL, TS, colour=eq)) +
  annotate("rect", xmin=2, xmax=135, ymin=min(ts_mid), ymax=max(ts_mid), fill="grey80", alpha=0.45) +
  annotate("text", x=2.1, y=max(ts_mid), label="survey TS range", hjust=0, vjust=-0.6, size=2.7, colour="grey40") +
  geom_line(linewidth=0.8) +
  geom_vline(xintercept=CATCH_MAX, linetype=2, colour="grey30") +
  annotate("text", x=CATCH_MAX, y=min(d4a$TS),
           label=sprintf("catch max\n%.1f cm",CATCH_MAX),
           hjust=1.1, vjust=0, size=2.5, colour="grey30") +
  scale_colour_manual(values=pal_tsl) + scale_x_log10() +
  labs(x="Total length (cm, log)", y="Target strength (dB)", colour=NULL, tag="a") +
  th + theme(legend.position="bottom")
d4b <- data.frame(eq=factor(tsl_short,levels=tsl_short), maxTL=maxTL)
p4b <- ggplot(d4b, aes(eq,maxTL,fill=eq)) + geom_col(width=0.7, show.legend=FALSE) +
  geom_hline(yintercept=CATCH_MAX, linetype=2, colour="grey30") +
  annotate("text", x=0.6, y=CATCH_MAX, label=sprintf("catch %.1f cm",CATCH_MAX),
           hjust=0, vjust=-0.5, size=2.4, colour="grey30") +
  scale_fill_manual(values=pal_tsl) + labs(x=NULL, y="Max inferred TL (cm)", tag="b") +
  th + theme(axis.text.x=element_text(angle=20,hjust=1))
d4c <- data.frame(eq=factor(tsl_short,levels=tsl_short), D=ks_D)
p4c <- ggplot(d4c, aes(eq,D,fill=eq)) + geom_col(width=0.7, show.legend=FALSE) +
  scale_fill_manual(values=pal_tsl) + labs(x=NULL, y="K-S D (acoustic vs catch)", tag="c") +
  th + theme(axis.text.x=element_text(angle=20,hjust=1))
lab4 <- sprintf("%s (D=%.2f)", tsl_short, ks_D)
d4d <- rbind(data.frame(len=gl_len, src="Catch (measured)"),
  do.call(rbind, lapply(seq_len(nrow(TSL)), function(i)
    data.frame(len=rep(TLen(TSL$m[i],TSL$b0[i],ts_mid), times=round(cls_tot)), src=lab4[i]))))
d4d$src <- factor(d4d$src, levels=c("Catch (measured)", lab4))
cols4d <- c("Catch (measured)"="black", setNames(as.character(pal_tsl[tsl_short]), lab4))
p4d <- ggplot(d4d, aes(len, colour=src)) + stat_ecdf(linewidth=0.7) +
  scale_x_log10() + scale_colour_manual(values=cols4d) +
  labs(x="Total length (cm, log)", y="Cumulative proportion", colour=NULL, tag="d") +
  th + theme(legend.position="inside", legend.position.inside=c(0.78,0.32), legend.text=element_text(size=6))
sav(p4a / (p4b|p4c) / p4d + plot_layout(heights=c(1.15,1,1.15)), "Fig4_TSlength_selection.png", 9, 11)

# =====================================================================
#  FIG 5 - sensitivity grid (Table 6), pooled biomass
# =====================================================================
grid <- outer(seq_len(nrow(TSL)), seq_len(nrow(LWR)),
              Vectorize(function(i,j) pooled_group_bio(TSL$m[i],TSL$b0[i],LWR$a[j],LWR$b[j])))
dimnames(grid) <- list(TSL$name, LWR$name)
gridout <- data.frame(TS_length_equation=TSL$name, frequency_kHz=TSL$freq, round(grid),
                      max_inferred_TL_cm=round(maxTL),
                      reproduces_catch_max=ifelse(maxTL>=CATCH_MAX,"yes","no"), check.names=FALSE)
write.csv(gridout, file.path(outdir,"Table6_sensitivity_grid.csv"), row.names=FALSE, fileEncoding="UTF-8")
ts_spread   <- max(grid[,"Community-fitted"])/min(grid[,"Community-fitted"])
lwr_spread_frou <- max(grid[sel,])/min(grid[sel,])
lwr_spread_max  <- max(apply(grid,1,function(r) max(r)/min(r)))
cat(sprintf("Grid: TS spread %.0f-fold ; LWR spread max %.2f-fold ; under Frouzova %.2f-fold\n",
            ts_spread, lwr_spread_max, lwr_spread_frou))

d5 <- data.frame(TSL=factor(rep(tsl_short,ncol(grid)),levels=tsl_short),
                 LWR=factor(rep(lwr_short,each=nrow(grid)),levels=lwr_short), t=as.vector(grid))
p5 <- ggplot(d5, aes(TSL,t,fill=LWR)) +
  geom_col(position=position_dodge(0.8), width=0.72, colour="grey25", linewidth=0.2) +
  scale_y_log10(labels=label_comma(), breaks=c(3e2,1e3,3e3,1e4,3e4,1e5,3e5)) +
  scale_fill_manual(values=pal_lwr) +
  labs(x="TS-length equation", y="Whole-lake biomass (t, log scale)", fill=NULL) +
  th + theme(legend.position="top")
sav(p5, "Fig5_sensitivity_grid.png", 8, 5)

# =====================================================================
#  FIG 6 - size-resolved biomass under Frouzova (Table 7), pooled method
# =====================================================================
glab <- BIO_LABELS
group_bio <- function(d, C, a, b){       # five pooled length-group totals (tonnes)
  pooled_group_bio(TSL$m[sel], TSL$b0[sel], a, b, d=d, C=C,
                   return_groups=TRUE)
}
pt <- list(ye=group_bio(den,Cnt,LWR$a[2],LWR$b[2]),
           comm=group_bio(den,Cnt,LWR$a[1],LWR$b[1]))
nb <- ceiling(nE/BLOCK); starts <- 1:(nE-BLOCK+1)
bootYe <- bootCo <- bootCR <- matrix(NA, BOOT_N, 5); bootTotYe <- bootTotCo <- numeric(BOOT_N)
set.seed(2025)                                         # isolate the spatial bootstrap
for (it in 1:BOOT_N){
  idx <- unlist(lapply(sample(starts, nb, replace=TRUE), function(s) s:(s+BLOCK-1)))[1:nE]
  # Recalculate both mean density and pooled TS composition in every replicate.
  ye <- group_bio(den[idx],Cnt[idx,,drop=FALSE],LWR$a[2],LWR$b[2])
  co <- group_bio(den[idx],Cnt[idx,,drop=FALSE],LWR$a[1],LWR$b[1])
  bootYe[it,]<-ye; bootCo[it,]<-co; bootCR[it,]<-(co-ye)/ye*100
  bootTotYe[it]<-sum(ye); bootTotCo[it]<-sum(co)
}
tYe<-sum(pt$ye); tCo<-sum(pt$comm)
ciYe<-quantile(bootTotYe,c(.025,.975)); ciCo<-quantile(bootTotCo,c(.025,.975))
crWL<-(tCo-tYe)/tYe*100; crWLci<-quantile((bootTotCo-bootTotYe)/bootTotYe*100,c(.025,.975))
CV<-sd(bootTotCo)/mean(bootTotCo)*100
qci<-function(M) apply(M,2,function(x) quantile(x,c(.025,.975),na.rm=TRUE))
ci_ye<-qci(bootYe); ci_co<-qci(bootCo); ci_cr<-qci(bootCR)
crGroup<-(pt$comm-pt$ye)/pt$ye*100
# Locate the CR zero crossing by interpolating between the pooled representative
# Frouzova lengths of adjacent fixed TS groups (descriptive output only).
all_counts <- colSums(Cnt, na.rm=TRUE)
group_rep_tl <- vapply(seq_along(glab), function(g){
  z <- TS_GROUP_ID==g
  mean_ts_g <- weighted.mean(ts_mid[z], all_counts[z], na.rm=TRUE)
  10^((mean_ts_g-TSL$b0[sel])/TSL$m[sel])
}, numeric(1))
cross_idx <- which(diff(sign(crGroup))!=0)
cross_group_cm <- if (length(cross_idx)) {
  i <- cross_idx[1]
  round(approx(crGroup[c(i,i+1)], group_rep_tl[c(i,i+1)], xout=0)$y, 1)
} else NA_real_

# Two-sided paired, centred-bootstrap P value.  The +1 correction prevents
# reporting an impossible P=0 with a finite number of bootstrap replicates.
paired_boot_p <- function(boot_diff, observed_diff){
  boot_diff <- boot_diff[is.finite(boot_diff)]
  if (!length(boot_diff) || !is.finite(observed_diff)) return(NA_real_)
  centred <- boot_diff-observed_diff
  (sum(abs(centred)>=abs(observed_diff))+1)/(length(centred)+1)
}
sig_symbol <- function(p) {
  ifelse(is.na(p), "NA",
    ifelse(p<0.001, "***", ifelse(p<0.01, "**", ifelse(p<0.05, "*", "n.s."))))
}
p_text <- function(p) {
  if (!is.finite(p)) return("P = NA")
  if (p<0.001) "P < 0.001" else sprintf("P = %.3f",p)
}
group_p <- vapply(seq_along(glab), function(g)
  paired_boot_p(bootCo[,g]-bootYe[,g], pt$comm[g]-pt$ye[g]), numeric(1))
whole_p <- paired_boot_p(bootTotCo-bootTotYe, tCo-tYe)
group_sig <- sig_symbol(group_p)
whole_sig <- sig_symbol(whole_p)

tab7 <- data.frame(length_group_cm=glab,
   Ye_borrowed_t=round(pt$ye,1), Ye_95CI=sprintf("%.0f to %.0f",ci_ye[1,],ci_ye[2,]),
   Community_t=round(pt$comm,1), Community_95CI=sprintf("%.0f to %.0f",ci_co[1,],ci_co[2,]),
   CR_percent=round(crGroup,1), CR_95CI=sprintf("%.1f to %.1f",ci_cr[1,],ci_cr[2,]),
   Bootstrap_P=round(group_p,4), Significance=group_sig,
   CI_excludes_zero=ifelse(ci_cr[1,]>0 | ci_cr[2,]<0,"yes","no"))
whole7 <- data.frame(
   length_group_cm="Whole lake",
   Ye_borrowed_t=round(tYe,1),
   Ye_95CI=sprintf("%.0f to %.0f",ciYe[1],ciYe[2]),
   Community_t=round(tCo,1),
   Community_95CI=sprintf("%.0f to %.0f",ciCo[1],ciCo[2]),
   CR_percent=round(crWL,1),
   CR_95CI=sprintf("%.1f to %.1f",crWLci[1],crWLci[2]),
   Bootstrap_P=round(whole_p,4),
   Significance=whole_sig,
   CI_excludes_zero=ifelse(crWLci[1]>0 | crWLci[2]<0,"yes","no"))
tab7 <- rbind(tab7, whole7)
write.csv(tab7, file.path(outdir,"Table7_size_resolved_biomass.csv"), row.names=FALSE, fileEncoding="UTF-8")

d6a <- rbind(data.frame(t=bootTotYe, m="Ye 2007"), data.frame(t=bootTotCo, m="Community-fitted"))
d6a$m <- factor(d6a$m, levels=c("Ye 2007","Community-fitted")); ytop <- max(quantile(bootTotYe,.99),quantile(bootTotCo,.99))
p6a <- ggplot(d6a, aes(m, t, fill=m)) +
  geom_violin(colour="grey20", linewidth=0.2, width=0.8) +
  stat_summary(fun=median, geom="point", shape=21, size=2.6, fill="white", colour="black") +
  annotate("segment", x=1, xend=2, y=ytop*1.03, yend=ytop*1.03, colour="grey30") +
  annotate("text", x=1.5, y=ytop*1.09, size=3.1,
           label=sprintf("CR = %+.1f%%; %s (%s)", crWL, p_text(whole_p), whole_sig)) +
  scale_fill_manual(values=c("Ye 2007"="#3C5488","Community-fitted"="#E64B35")) +
  scale_y_continuous(labels=label_comma()) +
  labs(x=NULL, y="Whole-lake biomass (t)", tag="a") + th + theme(legend.position="none")
d6b <- rbind(data.frame(grp=rep(glab,each=BOOT_N), t=as.vector(bootYe), m="Ye 2007"),
             data.frame(grp=rep(glab,each=BOOT_N), t=as.vector(bootCo), m="Community-fitted"))
d6b <- d6b[is.finite(d6b$t)&d6b$t>0,]; d6b$grp <- factor(d6b$grp, levels=glab)
d6b$m <- factor(d6b$m, levels=c("Ye 2007","Community-fitted"))
d6s <- data.frame(grp=factor(glab,levels=glab),
                  y=as.numeric(tapply(d6b$t,d6b$grp,function(x) quantile(x,.99))[glab])*1.6,
                  lab=group_sig)
p6b <- ggplot(d6b, aes(grp, t, fill=m)) +
  geom_violin(position=position_dodge(0.8), colour="grey30", linewidth=0.2, width=0.75, scale="width") +
  stat_summary(aes(group=m), fun=median, geom="point", position=position_dodge(0.8), shape=21, size=1.3, fill="white", colour="black") +
  geom_text(data=d6s, aes(grp,y,label=lab), inherit.aes=FALSE, size=3.6) +
  scale_fill_manual(values=c("Ye 2007"="#3C5488","Community-fitted"="#E64B35"), name=NULL) +
  scale_y_log10(labels=label_comma()) + labs(x="Length group (cm)", y="Biomass (t, log scale)", tag="b") +
  th + theme(legend.position="top")
d6c <- data.frame(grp=factor(glab,levels=glab), cr=crGroup, lo=ci_cr[1,], hi=ci_cr[2,])
d6c$sign <- ifelse(d6c$cr>=0,"Positive","Negative")
p6c <- ggplot(d6c, aes(grp, cr, colour=sign)) +
  geom_hline(yintercept=0, colour="grey50") +
  geom_segment(aes(xend=grp, y=0, yend=cr), linewidth=0.8) +
  geom_errorbar(aes(ymin=lo,ymax=hi), width=0.12, linewidth=0.4) + geom_point(size=3) +
  scale_colour_manual(values=c("Negative"="#3C5488","Positive"="#E64B35")) +
  labs(x="Length group (cm)", y="Correction rate (%)", tag="c") + th + theme(legend.position="none")
sav(p6a/p6b/p6c + plot_layout(heights=c(1,1.2,1)), "Fig6_size_resolved_biomass.png", 7, 11)

# =====================================================================
#  SUMMARY of every number reported in the manuscript
# =====================================================================
cat("\n=================  MANUSCRIPT NUMBERS  =================\n")
cat(sprintf("[3.1] dominance-size slope: numerical %.2f (95%% CI %.2f,%.2f, P=%.2f); biomass %.2f (95%% CI %.2f,%.2f, P=%.3f); cross %.0f cm\n",
    coef(snum)[2],confint(snum)[2,1],confint(snum)[2,2],summary(snum)$coef[2,4],
    coef(sbio)[2],confint(sbio)[2,1],confint(sbio)[2,2],summary(sbio)$coef[2,4], cross_cm))
cat(sprintf("[3.2] complete measured fish n=%d ; pooled LWR a=%.5f b=%.3f R2=%.3f ; species b range %.2f-%.2f ; Bayesian exponent %.2f (95%% CrI %.2f-%.2f)\n",
    nrow(cc), community_a, community_b, community_r2,
    min(res$b),max(res$b), comm_b, comm_ci[1], comm_ci[2]))
cat(sprintf("[3.2] species-slope ANCOVA F=%.2f, df=%d and %d, P=%s\n",
    anc_F, anc_df, anc_resid_df,
    ifelse(anc_P<0.001,"<0.001",sprintf("%.4f",anc_P))))
cat("[3.2] LWR accuracy vs measured catch (total bias %, median err %, log10-RMSE):\n"); print(accTab, row.names=FALSE)
cat(sprintf("[3.3] TS-equation spread %.0f-fold (%.0f to %.0f t) ; LWR spread <=%.2f-fold ; under Frouzova %.2f-fold\n",
    ts_spread, min(grid[,"Community-fitted"]), max(grid[,"Community-fitted"]), lwr_spread_max, lwr_spread_frou))
cat("      Table 6 grid (t):\n"); print(round(grid))
cat(sprintf("[3.4] max inferred TL: %s ; K-S D: %s\n",
    paste(round(maxTL),collapse="/"), paste(sprintf("%.2f",ks_D),collapse="/")))
cat(sprintf("[3.5] whole-lake Ye %.0f t (95%%CI %.0f-%.0f) ; Community %.0f t (95%%CI %.0f-%.0f)\n",
    tYe,ciYe[1],ciYe[2],tCo,ciCo[1],ciCo[2]))
cat(sprintf("[3.5] difference %+.1f%% (paired 95%%CI %+.1f to %+.1f, excludes 0: %s) ; paired bootstrap %s (%s) ; bootstrap CV %.1f%%\n",
    crWL,crWLci[1],crWLci[2], ifelse(crWLci[1]>0|crWLci[2]<0,"yes","no"),
    p_text(whole_p), whole_sig, CV))
cat(sprintf("[3.5] size-group CR from %+.1f%% to %+.1f%% ; crosses zero near %s cm\n",
    min(crGroup),max(crGroup), cross_group_cm))
cat("      Table 7:\n"); print(tab7, row.names=FALSE)
writeLines(capture.output(sessionInfo()), file.path(base, "session-info.txt"), useBytes=TRUE)
cat(sprintf("\nWrote tables to %s and figures to %s\n", outdir, figdir))
