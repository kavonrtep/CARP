#!/usr/bin/env Rscript
# make_workflow_figure.R -- the conceptual CARP schematic used in README.md.
#
#   scripts/make_workflow_figure.R [output_dir]      (default: figs)
#
# Produces figs/carp_workflow.{png,svg} from one drawing definition.
#
# This is deliberately NOT the same figure as figs/workflow_overview.svg, which
# scripts/make_workflow_diagram.py derives live from `snakemake --rulegraph` and
# which answers "which rules run in what order". This one answers "what does the
# pipeline do to a genome", and is hand-authored: it stays at the level of tools
# and concepts, never rule names or thresholds, because unlike the rulegraph it
# cannot be checked against the code and would otherwise drift.
#
# Visual vocabulary is shared with dante_ltr's README figures on purpose, so a
# reader moving between the two repositories is not relearning a language: a
# thin genome track spanning the full width, numbered stages down the page with
# the caption ABOVE its track (a left-hand caption gutter costs a fifth of the
# width and shrinks every glyph with it), earlier layers ghosted in later
# stages, and the Okabe-Ito palette with one rule -- a colour never means two
# things.
#
# Base R graphics only: R is already a runtime dependency (envs/tidecluster.yaml)
# so the figure needs nothing new to regenerate.

args   <- commandArgs(trailingOnly = TRUE)
OUTDIR <- if (length(args) >= 1) args[1] else "figs"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# --- palette (Okabe-Ito, softened) ----------------------------------------
# Fill is the hue at reduced strength, border the hue at full strength: flat
# blocks read as diagrams, a fill/border pair reads as an object.
COL <- list(
  ltr    = "#0072B2",  # blue            LTR retrotransposons  (DANTE_LTR)
  tir    = "#D55E00",  # vermillion      TIR DNA transposons   (DANTE_TIR)
  line   = "#CC79A7",  # reddish purple  LINEs                 (DANTE_LINE)
  tandem = "#009E73",  # bluish green    tandem repeats        (TideCluster)
  track  = "#D8D8D8",
  band   = "#F7F8FA",
  rule   = "#E3E6EA",
  ink    = "#14304A",
  mute   = "#7A8894",
  flag   = "#B22222",
  shadow = "#12325020"
)
ghost <- function(col, a = 0.14) adjustcolor(col, alpha.f = a)
soft  <- function(col, a = 0.80) adjustcolor(col, alpha.f = a)

AR <- 1.26          # y-units per x-unit for a visually round corner

# --- one genomic region, identical coordinates in every stage --------------
ELEM <- list(
  ltr  = list(c(11, 21), c(54, 64)),
  tir  = list(c(30, 36)),
  line = list(c(71, 77))
)
SAT   <- c(84, 97)
SAT_RM <- c(81.0, 97.0)          # RM extends the array beyond the clustered part
FRAG <- list(
  ltr  = list(c(4, 8), c(24, 28.5), c(41, 45), c(47.5, 52), c(66, 69.5)),
  tir  = list(c(38, 40), c(58.5, 60.5)),
  line = list(c(79, 81))
)
ORPHAN   <- list(c(45, 47.2))    # a domain with no element around it

# RepeatMasker also hits the ELEMENTS the library was built from -- the library
# consensi come from them, so they are the best matches in the genome. One hit
# runs past its element into flanking degraded sequence.
RM_HIT <- c(
  list(list(c(11, 25),   "ltr")),      # element + a degraded continuation
  lapply(ELEM$ltr[-1],  function(e) list(e, "ltr")),
  lapply(ELEM$tir,      function(e) list(e, "tir")),
  lapply(ELEM$line,     function(e) list(e, "line")),
  unlist(lapply(names(FRAG), function(k)
    lapply(FRAG[[k]], function(f) list(f, k))), recursive = FALSE)
)

# interval subtraction, so the reconciled track is provably non-overlapping
# rather than drawn to look that way
subtract <- function(iv, blockers) {
  out <- list(iv)
  for (b in blockers) {
    nxt <- list()
    for (p in out) {
      if (b[2] <= p[1] || b[1] >= p[2]) { nxt <- c(nxt, list(p)); next }
      if (b[1] > p[1]) nxt <- c(nxt, list(c(p[1], b[1])))
      if (b[2] < p[2]) nxt <- c(nxt, list(c(b[2], p[2])))
    }
    out <- nxt
  }
  Filter(function(p) p[2] - p[1] > 0.25, out)
}

# domains sit INSIDE the element body, clear of the terminal structure
domains_of <- function(e, kind, n) {
  pad <- if (kind == "ltr") 1.9 else if (kind == "tir") 1.9 else 0.4
  x0 <- e[1] + pad; x1 <- e[2] - pad
  w  <- (x1 - x0)/n
  lapply(seq_len(n) - 1, function(i) c(x0 + i*w + 0.25, x0 + (i+1)*w - 0.25))
}
DOM <- c(
  lapply(domains_of(ELEM$ltr[[1]],  "ltr",  5), function(d) list(d, "ltr")),
  lapply(domains_of(ELEM$ltr[[2]],  "ltr",  5), function(d) list(d, "ltr")),
  lapply(domains_of(ELEM$tir[[1]],  "tir",  1), function(d) list(d, "tir")),
  lapply(domains_of(ELEM$line[[1]], "line", 3), function(d) list(d, "line")),
  lapply(ORPHAN, function(d) list(d, "ltr"))
)

# --- primitives -----------------------------------------------------------
rrect <- function(x0, y0, x1, y1, r, col, border = NA, lwd = 1.4, shadow = TRUE, lty = 1) {
  rx <- min(r, (x1 - x0)/2)
  ry <- min(rx * AR, (y1 - y0)/2)
  k  <- 10
  arc <- function(cx, cy, t0, t1) {
    t <- seq(t0, t1, length.out = k)
    cbind(cx + rx*cos(t), cy + ry*sin(t))
  }
  P <- rbind(arc(x1-rx, y0+ry, -pi/2, 0),      # bottom-right
             arc(x1-rx, y1-ry,     0, pi/2),   # top-right
             arc(x0+rx, y1-ry,  pi/2, pi),     # top-left
             arc(x0+rx, y0+ry,    pi, 1.5*pi)) # bottom-left
  if (shadow) polygon(P[,1] + 0.12, P[,2] - 0.30, col = COL$shadow, border = NA)
  polygon(P[,1], P[,2], col = col, border = border, lwd = lwd, lty = lty)
}

track <- function(y) {
  rrect(0, y-0.42, 100, y+0.42, 0.42, COL$track, border = NA, shadow = FALSE)
}

# a complete element: soft body, saturated border, its own terminal structure
element <- function(x0, x1, y, h, kind, col, faded = FALSE) {
  if (faded) { rrect(x0, y-h, x1, y+h, 0.5, ghost(col), border = NA, shadow = FALSE); return(invisible()) }
  rrect(x0, y-h, x1, y+h, 0.5, soft(col, .82), border = col, lwd = 1.6)
  dk <- adjustcolor(col, red.f = .55, green.f = .55, blue.f = .55)
  if (kind == "ltr") {
    rrect(x0, y-h, x0+1.7, y+h, 0.5, dk, border = NA, shadow = FALSE)
    rrect(x1-1.7, y-h, x1, y+h, 0.5, dk, border = NA, shadow = FALSE)
  } else if (kind == "tir") {
    polygon(c(x0+0.1, x0+1.8, x0+1.8), c(y, y+h*0.94, y-h*0.94), col = dk, border = NA)
    polygon(c(x1-0.1, x1-1.8, x1-1.8), c(y, y+h*0.94, y-h*0.94), col = dk, border = NA)
  }
}

domain <- function(x0, x1, y, h, col, faded = FALSE) {
  tip <- min((x1-x0)*0.30, 1.1)
  if (faded) {
    polygon(c(x0, x1-tip, x1, x1-tip, x0), c(y-h, y-h, y, y+h, y+h),
            col = ghost(col, .30), border = NA)
    return(invisible())
  }
  polygon(c(x0+0.12, x1-tip+0.12, x1+0.12, x1-tip+0.12, x0+0.12),
          c(y-h-0.30, y-h-0.30, y-0.30, y+h-0.30, y+h-0.30), col = COL$shadow, border = NA)
  polygon(c(x0, x1-tip, x1, x1-tip, x0), c(y-h, y-h, y, y+h, y+h),
          col = soft(col, .55), border = col, lwd = 1.3)
}

satellite <- function(x0, x1, y, h, col, faded = FALSE, n = NULL) {
  if (is.null(n)) n <- max(4, round((x1-x0)/1.45))
  w <- (x1-x0)/n
  for (i in seq_len(n) - 1) {
    a <- x0 + i*w + 0.13; b <- x0 + (i+1)*w - 0.13
    if (faded) rrect(a, y-h, b, y+h, 0.28, ghost(col), border = NA, shadow = FALSE)
    else       rrect(a, y-h, b, y+h, 0.28, soft(col, .82), border = col, lwd = 1.1)
  }
}

# --- stage furniture ------------------------------------------------------
band <- function(ytop, ybot) rrect(-3.4, ybot, 103.4, ytop, 1.2, COL$band, border = NA, shadow = FALSE)

caption <- function(n, title, sub, y, right = NA, rcol = "#0F7B4F") {
  rrect(-2.6, y-1.45, 0.8, y+1.45, 1.7, COL$ink, border = NA, shadow = FALSE)
  text(-0.9, y, n, col = "white", cex = 0.72, font = 2)
  text(2.6, y, title, col = COL$ink, cex = 0.98, font = 2, adj = 0)
  if (!is.na(sub))
    text(2.6 + strwidth(title, cex = 0.98, font = 2) + 1.8, y, sub,
         col = COL$mute, cex = 0.70, adj = 0)
  if (!is.na(right)) {
    text(100, y, right, col = rcol, cex = 0.74, font = 2, adj = 1)
    text(100 - strwidth(right, cex = 0.74, font = 2) - 1.4, y, "✓", col = rcol, cex = 0.8, font = 2, adj = 1)
  }
  segments(2.6, y-2.3, 100, y-2.3, col = COL$rule, lwd = 1)
}

note <- function(x, y, txt, col = COL$mute, cex = 0.62, adj = 1) text(x, y, txt, col = col, cex = cex, adj = adj)

# a database cylinder -- the libraries are stores, not boxes; label goes outside.
# Fills are OPAQUE tints rather than alpha: the body, its bottom cap and the top
# cap overlap, and translucent fills composited into visible density bands.
tint <- function(col, p) {           # p = how far towards white
  v <- col2rgb(col)/255
  rgb(v[1]*(1-p) + p, v[2]*(1-p) + p, v[3]*(1-p) + p)
}
db_icon <- function(cx, cy, rx, h, col) {
  ry   <- rx * AR * 0.30
  body <- tint(col, 0.34)
  cap  <- tint(col, 0.58)
  seam <- tint(col, 0.62)
  arc  <- function(cyy, t0, t1, n = 48) {
    t <- seq(t0, t1, length.out = n)
    cbind(cx + rx*cos(t), cyy + ry*sin(t))
  }
  # one outline for body + bottom cap, so there is no seam where they meet
  B <- rbind(c(cx-rx, cy+h/2), arc(cy-h/2, pi, 2*pi), c(cx+rx, cy+h/2))
  polygon(B[,1] + 0.16, B[,2] - 0.38, col = COL$shadow, border = NA)
  polygon(B[,1], B[,2], col = body, border = col, lwd = 1.8)
  for (d in c(0.30, 0.62)) {          # stacked-disk seams, front half only
    A <- arc(cy + h/2 - d*h, pi, 2*pi, 40)
    lines(A[,1], A[,2], col = seam, lwd = 1.3)
  }
  T <- arc(cy + h/2, 0, 2*pi, 64)
  polygon(T[,1], T[,2], col = cap, border = col, lwd = 1.8)
}

# icons -------------------------------------------------------------------
doc_icon <- function(x, y, w, h, col, lab) {
  f <- w*0.30
  polygon(c(x+0.15, x+w-f+0.15, x+w+0.15, x+w+0.15, x+0.15),
          c(y+h-0.35, y+h-0.35, y+h-f-0.35, y-h-0.35, y-h-0.35), col = COL$shadow, border = NA)
  polygon(c(x, x+w-f, x+w, x+w, x), c(y+h, y+h, y+h-f, y-h, y-h),
          col = "white", border = col, lwd = 2)
  polygon(c(x+w-f, x+w-f, x+w), c(y+h, y+h-f, y+h-f), col = col, border = col)
  for (i in 1:3) segments(x+w*0.17, y+h*0.22-i*h*0.30, x+w*0.78, y+h*0.22-i*h*0.30, col = soft(col,.5), lwd = 1.5)
  text(x+w/2, y-h-2.3, lab, col = COL$ink, cex = 0.60, adj = 0.5)
}
pie_icon <- function(cx, cy, r, cols, fracs, lab) {
  a0 <- pi/2
  for (i in seq_along(fracs)) {
    a1 <- a0 - 2*pi*fracs[i]; aa <- seq(a0, a1, length.out = 44)
    polygon(c(cx, cx + r*cos(aa)), c(cy, cy + r*AR*sin(aa)),
            col = soft(cols[i], .88), border = "white", lwd = 1.4)
    a0 <- a1
  }
  text(cx, cy - r*AR - 2.3, lab, col = COL$ink, cex = 0.60, adj = 0.5)
}
wig_icon <- function(x, y, w, h, col, lab) {
  set.seed(7); xs <- seq(x, x+w, length.out = 36)
  ys <- y - h + abs(sin(seq(0, 6.2, length.out = 36)) * runif(36, .28, 1)) * 2*h
  for (i in seq_along(xs)) segments(xs[i], y-h, xs[i], ys[i], col = soft(col,.85), lwd = 2.2, lend = 1)
  segments(x, y-h, x+w, y-h, col = COL$track, lwd = 1.4)
  text(x+w/2, y-h-2.3, lab, col = COL$ink, cex = 0.60, adj = 0.5)
}

# --- the figure -----------------------------------------------------------
draw <- function() {
  par(mar = c(0.3, 0.3, 0.3, 0.3), xaxs = "i", yaxs = "i")
  plot(NA, xlim = c(-5, 105), ylim = c(-6, 103.5), axes = FALSE, xlab = "", ylab = "")
  H <- 1.7; RM <- 4.2

  # Stage baselines and band extents are explicit rather than a uniform pitch:
  # the stages are not the same height (4 has boxes, 6 has two lanes, 7 has
  # icons with labels beneath), so a single spacing rule puts a caption inside
  # the neighbouring band or a glyph through its own caption rule.
  Y <- c(s1 = 96, s2 = 83, s3 = 70, s4 = 51, s5 = 36, s6 = 18, s7 = 3.5)
  for (b in list(c(102.5, 88.0), c(76.5, 63.0), c(43.0, 27.5), c(10.5, -5.2)))
    band(b[1], b[2])

  ghosts <- function(y, sat = TRUE, dom = FALSE) {
    if (dom) for (d in DOM) domain(d[[1]][1], d[[1]][2], y, H, COL[[d[[2]]]], faded = TRUE)
    for (e in ELEM$ltr)  element(e[1], e[2], y, H, "ltr",  COL$ltr,  faded = TRUE)
    for (e in ELEM$tir)  element(e[1], e[2], y, H, "tir",  COL$tir,  faded = TRUE)
    for (e in ELEM$line) element(e[1], e[2], y, H, "line", COL$line, faded = TRUE)
    if (sat) satellite(SAT[1], SAT[2], y, H, COL$tandem, faded = TRUE)
  }

  ## 1 DANTE
  y <- Y["s1"]; caption(1, "DANTE", "protein domains — the evidence every structural call is built on", y + 4.6)
  track(y); for (d in DOM) domain(d[[1]][1], d[[1]][2], y, H, COL[[d[[2]]]])
  note(100, y - 4.4, "a domain with no element around it is itself an annotation layer")

  ## 2 complete elements
  y <- Y["s2"]; caption(2, "DANTE_LTR · DANTE_TIR · DANTE_LINE",
                        "complete elements, delimited by their own structure", y + 4.6)
  track(y)
  for (d in DOM) domain(d[[1]][1], d[[1]][2], y, H, COL[[d[[2]]]], faded = TRUE)
  for (e in ELEM$ltr)  element(e[1], e[2], y, H, "ltr",  COL$ltr)
  for (e in ELEM$tir)  element(e[1], e[2], y, H, "tir",  COL$tir)
  for (e in ELEM$line) element(e[1], e[2], y, H, "line", COL$line)
  lx <- 8; ly <- y - 4.6
  for (g in list(c("LTR-RT", COL$ltr), c("TIR", COL$tir), c("LINE", COL$line))) {
    rrect(lx, ly-0.85, lx+2.4, ly+0.85, 0.4, soft(g[2], .82), border = g[2], lwd = 1.3, shadow = FALSE)
    text(lx+3.2, ly, g[1], col = COL$mute, cex = 0.64, adj = 0)
    lx <- lx + 3.2 + strwidth(g[1], cex = 0.64) + 3.4
  }
  note(100, ly, "only a minority of the repeat content is intact enough to be found this way")

  ## 3 TideCluster
  y <- Y["s3"]; caption(3, "TideCluster", "tandem arrays, clustered into families", y + 4.6)
  track(y); ghosts(y, sat = FALSE, dom = TRUE)
  satellite(SAT[1], SAT[2], y, H, COL$tandem)
  note(mean(SAT), y - 4.4, "satellite array", col = COL$tandem, adj = 0.5)

  ## 4 libraries
  y <- Y["s4"]
  caption(4, "Repeat libraries", "built from what was just found — the one step that leaves the genome", y + 9.5)
  db_icon(16, y, 5.2, 4.6, COL$ltr)
  text(23.5, y + 1.2, "dispersed repeat library", col = COL$ink, cex = 0.72, font = 2, adj = 0)
  text(23.5, y - 0.9, "consensi of the complete LTR / TIR / LINE elements",
       col = COL$mute, cex = 0.60, adj = 0)
  db_icon(90, y, 5.2, 4.6, COL$tandem)
  text(82.5, y + 1.2, "tandem dimer library", col = COL$ink, cex = 0.72, font = 2, adj = 1)
  text(82.5, y - 0.9, "consensus monomers, doubled", col = COL$mute, cex = 0.60, adj = 1)
  # arrows start under the layers that feed each library and reach into stage 5
  arrows(16, y + 6.9, 16, y + 4.5, length = 0.07, col = COL$ltr,    lwd = 1.9)
  arrows(90, y + 6.9, 90, y + 4.5, length = 0.07, col = COL$tandem, lwd = 1.9)
  arrows(16, y - 4.6, 16, y - 7.6, length = 0.07, col = COL$ltr,    lwd = 1.9)
  arrows(90, y - 4.6, 90, y - 7.6, length = 0.07, col = COL$tandem, lwd = 1.9)

  ## 5 RepeatMasker
  y <- Y["s5"]; caption(5, "RepeatMasker", "searches both libraries back against the whole genome", y + 4.6)
  track(y); ghosts(y)
  for (h in RM_HIT)
    rrect(h[[1]][1], y-RM-H*0.8, h[[1]][2], y-RM+H*0.8, 0.4,
          soft(COL[[h[[2]]]], .80), border = COL[[h[[2]]]], lwd = 1.2)
  satellite(SAT_RM[1], SAT_RM[2], y-RM, H*0.8, COL$tandem)
  note(100, y - RM - 3.6,
       "the source elements match too — plus the degraded copies, which are most of the repeat content")

  ## 6 reconciliation -- two lanes: what came from structure, what from similarity
  y <- Y["s6"]
  caption(6, "Reconciliation", "layers overlap; tier priority resolves every base to one call",
          y + 7.2, right = "one non-overlapping annotation")
  track(y)
  ys <- y + 2.8; yr <- y - 2.8
  blockers <- c(ELEM$ltr, ELEM$tir, ELEM$line, list(SAT), ORPHAN)

  # upper lane: structural evidence wins wherever it exists
  for (e in ELEM$ltr)  element(e[1], e[2], ys, H, "ltr",  COL$ltr)
  for (e in ELEM$tir)  element(e[1], e[2], ys, H, "tir",  COL$tir)
  for (e in ELEM$line) element(e[1], e[2], ys, H, "line", COL$line)
  satellite(SAT[1], SAT[2], ys, H, COL$tandem)
  for (o in ORPHAN) domain(o[1], o[2], ys, H, COL$ltr)   # a domain is a call of its own

  # lower lane: similarity hits, minus everything the upper lane already claimed
  for (h in RM_HIT) for (p in subtract(h[[1]], blockers))
    rrect(p[1], yr-H*0.85, p[2], yr+H*0.85, 0.4,
          soft(COL[[h[[2]]]], .80), border = COL[[h[[2]]]], lwd = 1.2)
  for (p in subtract(SAT_RM, blockers))
    satellite(p[1], p[2], yr, H*0.85, COL$tandem)

  note(100, y - 6.4, "upper: structure-based calls    ·    lower: similarity-based calls")

  ## 7 outputs
  y <- Y["s7"]; caption(7, "Outputs", NA, y + 6.8)
  doc_icon(11, y, 8.0, 3.7, COL$ink,  "GFF3 annotation")
  wig_icon(32, y + 3.7, 14, 3.7, COL$ltr, "density BigWigs")
  pie_icon(62, y, 4.1, c(COL$ltr, COL$tir, COL$line, COL$tandem, COL$mute),
           c(.42, .14, .10, .16, .18), "summary statistics")
  doc_icon(80, y, 8.0, 3.7, COL$flag, "interactive HTML report")
}

render <- function(path, dev, ...) {
  dev(path, ...); draw(); invisible(dev.off()); message("wrote ", path)
}
render(file.path(OUTDIR, "carp_workflow.png"), png,
       width = 1800, height = 1540, res = 150, bg = "white")
render(file.path(OUTDIR, "carp_workflow.svg"), svg,
       width = 12, height = 10.3, bg = "white")
