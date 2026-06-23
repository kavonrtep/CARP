#!/usr/bin/env Rscript
main_dir <- commandArgs(TRUE)[1]
output_pdf <- commandArgs(TRUE)[2]
# for testing:
if (FALSE){
  main_dir <- "/mnt/ceph/454_data/DToL/Henderson_paper/tests/container_repeat_annot/JunEffu_240619"
}
suppressPackageStartupMessages({
  library(rtracklayer)
})

# Robustness: on ANY uncaught error (e.g. data too large to render, an
# allocation failure, a malformed input) emit a one-page placeholder PDF and
# exit 0, so the Snakemake rule still gets a valid output instead of failing or
# leaving a truncated file. Per-panel failures are handled separately below
# (placeholder panel); this global handler catches everything else.
write_placeholder_pdf <- function(out, msg){
  try(grDevices::graphics.off(), silent = TRUE)  # close any half-open device
  try({
    pdf(out, width = 8, height = 6)
    plot.new()
    text(0.5, 0.55, "Summary plots not rendered", cex = 1.6, col = "grey25")
    text(0.5, 0.42, msg, cex = 0.8, col = "grey55")
    dev.off()
  }, silent = TRUE)
}
options(error = function(){
  write_placeholder_pdf(output_pdf, geterrmessage())
  quit(save = "no", status = 0)
})

smooth_score <- function(x, N_for_mean = 10){
  # extend the score in each direction by N_for_mean-1 zeros
  sc <- c(rep(0, N_for_mean-1), x, rep(0, N_for_mean-1))
  sc_smooth <- filter(sc, rep(1/N_for_mean, N_for_mean), sides=2)
  # remove the first N_for_mean-1 and the last N_for_mean-1 elements
  x <- sc_smooth[(N_for_mean):(length(sc_smooth)-N_for_mean+1)]
  x[1:N_for_mean-1] <- NA
  x[(length(x)-N_for_mean+2):length(x)] <- NA
  x
}

plot_tracks <- function(main_tracks, SL,
                        loess_50_threshold = 50, loess_100_threshold = 100,reverse = FALSE,
                        col1 = "#00000050", col2 = "#FF000050", col3 = "#FF0000FF",
                        empty_msg = "no tracks available") {
  # Drop NULL/empty tracks and bail out gracefully on an empty set, so a
  # genome lacking a whole category (e.g. no LTR lineages or no satellites)
  # yields a labelled placeholder page instead of crashing mid-PDF.
  main_tracks <- main_tracks[!vapply(main_tracks, is.null, logical(1))]
  main_tracks <- main_tracks[vapply(main_tracks, function(x) length(x) > 0, logical(1))]
  if (length(main_tracks) == 0){
    plot.new(); text(0.5, 0.5, empty_msg, cex = 1.2, col = "grey40")
    return(invisible(NULL))
  }
  ymax <- max(sapply(main_tracks, function(x) max(x$score)), na.rm = TRUE)
  if (!is.finite(ymax) || ymax <= 0) ymax <- 1
  par(mar=c(5,1,1,1))
  if (reverse){
    ylims <- c(length(main_tracks)*1.1 * ymax, 0)
  }else{
    ylims <- c(0, length(main_tracks)*1.1 * ymax)
  }


  plot(0,0, type='n', xlim=c(0,sum(SL)*1.2), ylim=ylims,
       xlab="Genome position", ylab="", axes=FALSE)

  main_tracks_df_list <- list()
  for (i in seq_along(main_tracks)){
    x <- as.data.frame(main_tracks[[i]])
    x$track <- names(main_tracks)[i]
    x$position <- (x$start + x$end) / 2
    main_tracks_df_list[[i]] <- x
  }
  names(main_tracks_df_list) <- names(main_tracks)
  SN_offsets <- cumsum(c(0, SL))
  names(SN_offsets) <- names(SL)
  yoffest <- 0

  for (tt in names(main_tracks)){
    xx <- main_tracks_df_list[[tt]]
    for (sn in names(SL)){
      x_part <- xx[xx$seqnames == sn,]
      xcoords <- x_part$position + SN_offsets[sn]
      ycoords <- x_part$score + yoffest

      # smooth line using loess
      if (length(xcoords) > loess_50_threshold){
        ycoords_smooth <- smooth_score(ycoords, 50)
        points(xcoords, ycoords, col = col1, pch=18, cex=0.2, type = "l")
        points(xcoords, ycoords_smooth, col = col2, pch=18, cex=0.2, type = "l")
      }
      if (length(xcoords) > loess_100_threshold){
        ycoords_smooth2 <- smooth_score(ycoords, 100)
        points(xcoords, ycoords_smooth2, col = col3, pch=18, cex=0.2, type = "l")
      }
    }
    # add labels on side - track names
    text(max(SN_offsets), yoffest + 0, paste(tt,"\n "), pos = 4)
    yoffest <- yoffest + 1.1 * ymax
  }

  axis(1, at = (head(SN_offsets, -1) + tail(SN_offsets, -1))/2, labels = names(SL),
       # not ticks and lines
       lwd = 0, cex.axis = 0.5, line = -1, las = 2, adj = 1)

  # use segments instead of abline to draw horizontal lines
  for (i in 0:length(main_tracks)){
    segments(0, i * 1.1 * ymax , sum(SL), i*1.1 * ymax, col = "grey")
  }

  # vertical lines
  # no lines if length os SL is 1
  if (length(SL) > 1){
    for (i in 1:(length(SL)-1)){
      segments(SN_offsets[i], 0, SN_offsets[i], length(main_tracks)*1.1, col = "grey")
    }
  }
}



SL <- readRDS(paste0(main_dir, "/genome_seqlengths.rds"))
dir_100k_RA <- paste0(main_dir,"/Repeat_density_by_class_bigwig/100k")
# sort by lengths
SL <- SL[order(SL, decreasing = TRUE)]
widths <- SL/sum(SL)*100

list_of_tracks <-  c(
  Mobile_elements = paste0(dir_100k_RA, "/Mobile_elements_100k.bw"),
  Low_complexity = paste0(dir_100k_RA, "/Low_complexity_100k.bw"),
  Simple_repeats = paste0(dir_100k_RA, "/Simple_repeat_100k.bw"),
  rDNA = paste0(dir_100k_RA, "/rDNA_100k.bw"),
  Ty1_Copia = paste0(dir_100k_RA, "/All_Ty1_Copia_100k.bw"),
  Ty3_Gypsy = paste0(dir_100k_RA, "/All_Ty3_Gypsy_100k.bw"),
  Tandem_repeats_TC = paste0(main_dir, "/Tandem_repeats_TideCluster_100k.bw")
)
# load only existing files
main_tracks <- list()
for (name in names(list_of_tracks)){
  if (file.exists(list_of_tracks[[name]])){
    main_tracks[[name]] <- import(list_of_tracks[[name]])
  }
}


# conver to data.frame for ggplot, uses seqnames, mean(start+end), score and track name
main_tracks_df_list <- list()
for (i in seq_along(main_tracks)){
  x <- as.data.frame(main_tracks[[i]])
  x$track <- names(main_tracks)[i]
  x$position <- (x$start + x$end) / 2
  main_tracks_df_list[[i]] <- x
}
names(main_tracks_df_list) <- names(main_tracks)


# plot lineges
RM_dir <- paste0(main_dir, "/Repeat_density_by_class_bigwig/100k")
lineages_file <- dir(RM_dir, pattern = "LTR.Ty")
# keep only defined lineages
lineages <- c('Ale', 'Alesia', 'Angela', 'Bianca', 'Bryco', 'Gymco-I', 'Gymco-II',
              'Gymco-III', 'Gymco-IV', 'Ikeros', 'Ivana', 'Lyco', 'Osser', 'SIRE',
              'TAR', 'Tork', 'Chlamyvir', 'chromo-unclass', 'CRM', 'Galadriel',
              'Reina', 'Tcn1', 'Tekay', 'Athila', 'Ogre', 'Retand', 'TatI',
              'TatII', 'TatIII', 'Phygy', 'Selgy')
# keep only files with defined lineages string
lineages_file <- lineages_file[grep(paste0(lineages, collapse = "|"), lineages_file)]

lineages_bw <- lapply(lineages_file, function(x)import(paste0(RM_dir, "/", x)))
names(lineages_bw) <- gsub(".+[.]", "",
                           gsub("_100.+", "", lineages_file))





# plot major satellites

read_csv_or_empty <- function(path){
  if (file.exists(path)) {
    tryCatch(read.table(path, header = TRUE, sep = "\t"),
             error = function(e) data.frame())
  } else {
    data.frame()
  }
}
monomer_estimates_kite <- read_csv_or_empty(
  paste0(main_dir, "/TideCluster/default/TideCluster_kite/monomer_size_top3_estimats.csv"))

monomer_best_estimage <- read_csv_or_empty(
  paste0(main_dir, "/TideCluster/default/TideCluster_kite/monomer_size_best_estimate_stat.csv"))

major_bigwig <- dir(
  paste0(main_dir,
         "/Tandem_repeats_TideCluster_split_by_family_bigwig/100k"),
  pattern = "bw")

trc_name <- gsub("_100k.bw", "", major_bigwig)
trc_index <- as.numeric(gsub("TRC_", "", trc_name))

# reorder by index
major_bigwig <- major_bigwig[order(trc_index)]
trc_name <- trc_name[order(trc_index)]
trc_index <- trc_index[order(trc_index)]

# read first max 20 satellites (seq_len handles the no-satellite case: N == 0)
N <- min(20, length(major_bigwig))
trc_bw <- list()
have_monomer <- all(c("position", "TRC_ID") %in% names(monomer_best_estimage))
for (i in seq_len(N)){
  label <- trc_name[i]
  monomer_size <- NA
  if (have_monomer){
    monomer_size <- names(sort(table(monomer_best_estimage$position[monomer_best_estimage$TRC_ID == label]),
                          decreasing = TRUE)[1])
  }
  lab <- if (is.null(monomer_size) || is.na(monomer_size)) label else paste0(label, " (", monomer_size, "bp)")
  trc_bw[[lab]] <- import(paste0(main_dir, "/Tandem_repeats_TideCluster_split_by_family_bigwig/100k/",
                                     major_bigwig[[i]]))
}




pdf(output_pdf, width = 14, height = 7)
on.exit(dev.off(), add = TRUE)   # always finalise a valid PDF
# Each panel is isolated: an error in one (e.g. an empty category) draws a
# placeholder and does not abort the others or truncate the file.
safe_panel <- function(expr, msg){
  tryCatch(expr, error = function(e){
    message("summary_plots: ", msg, " panel failed: ", conditionMessage(e))
    plot.new(); text(0.5, 0.5, paste0(msg, " unavailable"), cex = 1.2, col = "grey40")
  })
}
safe_panel(plot_tracks(main_tracks, SL, empty_msg = "no repeat-class tracks"),
           "repeat-class")
safe_panel(plot_tracks(lineages_bw, SL, empty_msg = "no LTR lineage tracks"),
           "LTR lineage")
safe_panel(plot_tracks(rev(trc_bw), SL, empty_msg = "no satellite (TRC) tracks"),
           "satellite")


