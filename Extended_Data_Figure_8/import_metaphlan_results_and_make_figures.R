library(phyloseq)
library(microbiome)
library(tidyverse)
library(vegan)
library(ggpubr)
library(patchwork)

outputdir <- "./tmp/"

# --- Helper functions ---

aggregate_topn <- function(x, level, detection=0, prevalence=0.5, topn=10, toi=NULL) {
  x <- aggregate_taxa(x, level)
  tax <- phyloseq::tax_table(x)
  top_prev <- names(microbiome::prevalence(x, detection, sort=TRUE))[1:topn]
  top_prev <- unique(c(top_prev[top_prev %!in% toi][1:(topn - length(toi))], toi))
  inds <- which(!(rownames(tax) %in% top_prev))
  tax[inds, level] <- "other"
  phyloseq::tax_table(x) <- tax
  tt <- phyloseq::tax_table(x)[, level]
  phyloseq::tax_table(x) <- phyloseq::tax_table(tt)
  aggregate_taxa(x, level)
}

myotudf <- function(ps, rankname="Family", vars_of_interest=NULL) {
  ps <- transform_sample_counts(ps, function(OTU) OTU / sum(OTU))
  samdf <- data.frame(phyloseq::sample_data(ps)) %>%
    tibble::rownames_to_column(var="mysuperuniquesampleID") %>%
    mutate(mysuperuniquesampleID=make.names(mysuperuniquesampleID))
  otudf <- data.frame(otu_table(ps), stringsAsFactors=FALSE)
  taxdf <- data.frame(phyloseq::tax_table(ps), stringsAsFactors=FALSE)
  df <- merge(taxdf, otudf, by=0)
  allranks <- rank_names(ps)[1:which(rank_names(ps) == rankname)]
  df %>%
    group_by_at(allranks) %>%
    summarise_if(is.numeric, .funs=(sum)) %>%
    ungroup() %>%
    gather(mysuperuniquesampleID, val, (length(allranks) + 1):ncol(.)) %>%
    pivot_wider(id_cols="mysuperuniquesampleID", names_from=all_of(allranks), values_from="val", names_sep=";") %>%
    inner_join(., samdf %>% dplyr::select(mysuperuniquesampleID, all_of(vars_of_interest)), by="mysuperuniquesampleID") %>%
    dplyr::rename(SID=mysuperuniquesampleID) %>%
    dplyr::select(SID, all_of(vars_of_interest), everything())
}

myadonis2 <- function(ps=NULL, pretext=NULL, adoby="terms", test_groups=c("Group"), strata=NULL,
                      outputdir="pcoa/", dist="bray", gdist=NULL, name=NULL, norandom=FALSE, nofile=FALSE) {
  if (!file.exists(outputdir) & !nofile) system(paste("mkdir", outputdir))
  phyloseq::sample_data(ps)$Group <- unlist(phyloseq::sample_data(ps)[, test_groups[1]])
  ps <- phyloseq::subset_samples(ps, !is.na(Group))
  ps <- transform_sample_counts(ps, function(OTU) OTU / sum(OTU))
  if (is.null(gdist)) {
    if (dist == "gunifrac") {
      bc.ps <- gUniFrac(ps)
      save(bc.ps, file=paste0(outputdir, dist))
    } else if (dist == "wunifrac") {
      tree <- phy_tree(ps)
      OTU <- as.matrix(otu_table(ps))
      if (taxa_are_rows(ps)) OTU <- t(OTU)
      unifracs <- GUniFrac::GUniFrac(OTU, tree=tree)$unifracs
      bc.ps <- as.dist(unifracs[, , "d_1"])
      save(bc.ps, file=paste0(outputdir, dist))
    } else if (dist == "gunifrac_vaw") {
      tree <- phy_tree(ps)
      OTU <- as.matrix(otu_table(ps))
      if (taxa_are_rows(ps)) OTU <- t(OTU)
      unifracs <- GUniFrac::GUniFrac(OTU, tree=tree)$unifracs
      bc.ps <- as.dist(unifracs[, , "d_VAW"])
      save(bc.ps, file=paste0(outputdir, dist))
    } else {
      bc.ps <- phyloseq::distance(ps, method=dist)
    }
  } else {
    bc.ps <- gdist
  }
  df <- data.frame(phyloseq::sample_data(ps))
  if (!norandom) df$random1 <- sample(letters[1:2], nsamples(ps), TRUE)
  form1 <- if (!norandom) {
    formula(paste0("bc.ps ~ ", paste0(test_groups, collapse="+"), "+random1"))
  } else {
    formula(paste0("bc.ps ~ ", paste0(test_groups, collapse="+")))
  }
  if (is.null(strata)) {
    result1 <- adonis2(form1, data=df, by=adoby)
  } else {
    perm <- how(nperm=999, 2)
    setBlocks(perm) <- with(df, strata)
    result1 <- adonis2(form1, data=df, permutations=perm)
  }
  if (!nofile) {
    sink(paste0(outputdir, pretext, name, paste0(test_groups, collapse="&"), "_", dist, "_adonis2.txt"))
    print(form1)
    print(result1)
    sink()
    sink(paste0(outputdir, pretext, name, paste0(test_groups, collapse="&"), "_", dist, "_anosim.txt"))
    print(test_groups)
    for (var in test_groups) {
      if (is.null(strata)) {
        text <- capture.output(anosim(bc.ps, df[, var]))
      } else {
        text <- capture.output(anosim(bc.ps, df[, var], strata=df[, strata]))
      }
      text[3] <- gsub("..., var.", var, text[3])
      print(cat(text, sep="\n"), na.print="")
    }
    sink()
  }
  return(list(bc.ps, df, form1, result1))
}

mypcoa <- function(ps=NULL, plot_type="NULL", colors=NULL, legtex="", shapetex="", name=NULL,
                   shape=NULL, shape_keep_group=NULL, label=NULL, lpos="bottom", lbox="horizontal",
                   pretext=NULL, test_group="Group", nocap=FALSE, outputdir="pcoa/", dist="bray",
                   method="NMDS", type="both", distname=NULL, gdist=NULL,
                   norandom=TRUE, cex=1) {
  if (!file.exists(outputdir)) system(paste("mkdir", outputdir))
  sample_names(ps) <- gsub("^(\\d+)", "X\\1", sample_names(ps))
  if (is.null(distname)) {
    distname <- switch(dist,
      bray     = paste0("Bray-Curtis ", method),
      wunifrac = paste0("Weighted UniFrac ", method),
      gunifrac = paste0("Generalized UniFrac ", method),
      unifrac  = paste0("UniFrac ", method),
      jsd      = paste0("Jensen-Shannon Divergence ", method),
      ""
    )
  }
  phyloseq::sample_data(ps)$Group <- unlist(phyloseq::sample_data(ps)[, test_group])
  ps <- phyloseq::subset_samples(ps, !is.na(Group))
  ps <- transform_sample_counts(ps, function(OTU) OTU / sum(OTU))
  if (is.null(gdist)) {
    if (dist %in% c("gunifrac", "wunifrac", "gunifrac_vaw")) {
      tree <- phy_tree(ps)
      OTU <- as.matrix(otu_table(ps))
      if (taxa_are_rows(ps)) OTU <- t(OTU)
      unifracs <- GUniFrac::GUniFrac(OTU, tree=tree)$unifracs
      bc.ps <- as.dist(unifracs[, , switch(dist, gunifrac="d_0.5", wunifrac="d_1", gunifrac_vaw="d_VAW")])
      saveRDS(bc.ps, file=paste0(outputdir, pretext, dist, ".rds"))
    } else {
      bc.ps <- phyloseq::distance(ps, method=dist)
    }
  } else {
    bc.ps <- gdist
  }
  if (method == "CAP") {
    ord.NMDS.bc <- ordinate(ps, method="CAP", distance=bc.ps, formula=~Group)
  } else {
    ord.NMDS.bc <- ordinate(ps, method=method, distance=bc.ps)
  }
  if (is.null(shape_keep_group)) {
    p <- plot_ordination(ps, ord.NMDS.bc, color="Group", shape=shape, title="") + geom_point(size=cex)
  } else {
    phyloseq::sample_data(ps)$shape2 <- unlist(phyloseq::sample_data(ps)[, shape_keep_group])
    p <- plot_ordination(ps, ord.NMDS.bc, color="Group", title="") + geom_point(size=cex * 1.5, aes(shape=shape2))
    p$layers[[1]]$aes_params$size <- 0
  }
  if (!is.null(label)) p <- p + geom_text(aes_string(label=label), size=cex * 2, vjust=1.5)
  if (plot_type == "classic") p <- p + theme_classic2()
  if (plot_type == "bw")      p <- p + theme_bw()
  if (type != "none") {
    p <- p + switch(type,
      norm         = list(stat_ellipse(level=0.95, type="norm", aes(color=Group, fill=Group), alpha=0.1, geom="polygon"), stat_stars(alpha=0.3)),
      norm_no_star = list(stat_ellipse(level=0.95, type="norm", aes(color=Group, fill=Group), alpha=0.1, geom="polygon")),
      t            = list(stat_ellipse(level=0.95, type="t",    aes(color=Group, fill=Group), alpha=0.1, geom="polygon"), stat_stars(alpha=0.3)),
      t_no_star    = list(stat_ellipse(level=0.95, type="t",    aes(color=Group, fill=Group), alpha=0.1, geom="polygon")),
      both         = list(stat_ellipse(type="norm", linetype=2), stat_stars(alpha=0.3), stat_conf_ellipse(level=0.95, aes(color=Group, fill=Group), alpha=0.1, geom="polygon")),
      both_no_star = list(stat_ellipse(type="norm", linetype=2), stat_conf_ellipse(level=0.95, aes(color=Group, fill=Group), alpha=0.1, geom="polygon")),
      conf         = list(stat_conf_ellipse(level=0.95, aes(shape=Group, fill=Group), alpha=0.1, geom="polygon")),
      confstar     = list(stat_conf_ellipse(level=0.95, aes(color=Group, fill=Group), alpha=0.1, geom="polygon"), stat_stars(alpha=0.3)),
      star         = list(stat_stars(alpha=0.6)),
      NULL
    )
    caption <- switch(type,
      norm         = ,
      norm_no_star = paste0(distname, "\n(95% confidence levels assuming normal (\u2013\u2013) distribution)"),
      t            = ,
      t_no_star    = paste0(distname, "\n(95% confidence levels assuming multivariate t- (\u2013\u2013) distribution)"),
      both         = ,
      both_no_star = paste0(distname, "\n(95% confidence levels assuming normal (---) distribution and 95% confidence ellipses (\u2013\u2013))"),
      conf         = ,
      confstar     = paste0(distname, "\n(95% confidence ellipses (\u2013\u2013))"),
      NULL
    )
    if (!is.null(caption)) p <- p + labs(caption=caption)
  }
  if (!is.null(colors)) p <- p + scale_color_manual(values=colors) + scale_fill_manual(values=colors)
  p <- p +
    theme(legend.position=lpos, legend.box=lbox, plot.caption=element_text(size=rel(0.5))) +
    labs(fill=legtex, color=legtex, shape=shapetex)
  if (nocap) p <- p + labs(caption=NULL)
  result1 <- myadonis2(ps=ps, test_groups=c("Group", shape), dist=dist, gdist=bc.ps, norandom=norandom, nofile=TRUE)[[4]]
  form1 <- if (!norandom) {
    formula(paste0("bc.ps ~ ", paste0(c(test_group, shape), collapse="+"), "+random1"))
  } else {
    formula(paste0("bc.ps ~ ", paste0(c(test_group, shape), collapse="+")))
  }
  p %>% ggsave(plot=., filename=paste0(outputdir, pretext, name, test_group, "_", dist, "_", method, ".pdf"),
               width=297 / sqrt(2), height=210 / sqrt(2), units="mm")
  sink(paste0(outputdir, pretext, name, test_group, "_", dist, "_stat.txt"))
  print(form1)
  print(result1)
  cat("\n")
  print(data.frame(phyloseq::sample_data(ps), stringsAsFactors=FALSE) %>% count(Group))
  sink()
  return(list(p))
}

metaphlanToPhyloseq <- function(
    tax,
    metadat=NULL,
    simplenames=TRUE,
    roundtointeger=FALSE,
    split="|"){
  ## tax is a matrix or data.frame with the table of taxonomic abundances, rows are taxa, columns are samples
  ## metadat is an optional data.frame of specimen metadata, rows are samples, columns are variables
  ## if simplenames=TRUE, use only the most detailed level of taxa names in the final object
  ## if roundtointeger=TRUE, values will be rounded to the nearest integer
  xnames = rownames(tax)
  shortnames = gsub(paste0(".+\\", split), "", xnames)
  if(simplenames){
    rownames(tax) = shortnames
  }
  if(roundtointeger){
    tax = round(tax * 1e4)
  }
  x2 = strsplit(xnames, split=split, fixed=TRUE)
  taxmat = matrix(NA, ncol=max(sapply(x2, length)), nrow=length(x2))
  colnames(taxmat) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")[1:ncol(taxmat)]
  rownames(taxmat) = rownames(tax)
  for (i in 1:nrow(taxmat)){
    taxmat[i, 1:length(x2[[i]])] <- x2[[i]]
  }
  taxmat = gsub("[a-z]__", "", taxmat)
  taxmat = phyloseq::tax_table(taxmat)
  otutab = phyloseq::otu_table(tax, taxa_are_rows=TRUE)
  if(is.null(metadat)){
    res = phyloseq::phyloseq(taxmat, otutab)
  }else{
    res = phyloseq::phyloseq(taxmat, otutab, phyloseq::sample_data(metadat))
  }
  return(res)
}

# --- Import ---

read_metaphlan <- function(path, id_pattern) {
  abundance <- readr::read_delim(path, delim = "\t", trim_ws = TRUE, skip = 1) %>%
    tibble::column_to_rownames(colnames(.)[1])
  ps <- metaphlanToPhyloseq(abundance)
  sample_names(ps) <- gsub(id_pattern, "Unknown_\\1", sample_names(ps))
  ps
}

phyloseqin <- merge_phyloseq(
  read_metaphlan("./data/merged_abundance_table1.tsv", ".*_(BM\\d{3}-....-\\d{6})_.*"),
  read_metaphlan("./data/merged_abundance_table2.tsv", ".*_(CH\\d{3}-001R\\d{4})_.*")
)

taxa_names(phyloseqin) <- gsub("_group", "", taxa_names(phyloseqin))
ps <- prune_taxa(phyloseqin, taxa = grep("t__", taxa_names(phyloseqin), value = TRUE))

day_levels <- c("Donor", "Baseline", "Day 8", "Day 15-19", "Day 29-43")

map <- readr::read_delim("./data/attributes_table.tsv") %>%
  left_join(readxl::read_xlsx("./data/seq_days.xlsx"), by = join_by("sid" == "SeqID")) %>%
  transmute(
    SeqID            = sid,
    Pat_ID           = recode(Patient, "1" = "Pat. 1", "2" = "Pat. 2", "3" = "Pat. 3"),
    sample_type_spec = case_when(
      env_local_scale == "ENVO:00002003" ~ "stool",
      env_local_scale == "ENVO:00002047" ~ "urine",
      env_local_scale == "ENVO:02000047" ~ "vaginal swab"
    ),
    day              = day,
    day_cat          = factor(gsub("Day", "Day ", timepoint), levels = day_levels)
  )

sample_data(ps) <- map %>% tibble::column_to_rownames("SeqID")

# --- Figures ---

patient_colors <- c("Pat. 1" = "#2FB3CA", "Pat. 2" = "#B894C0", "Pat. 3" = "#008000")
focal_patients <- names(patient_colors)
sample_type_labels <- c(stool = "Stool", urine = "Urine", "vaginal swab" = "Vswab")

# Fig. 9A: E. coli relative abundance across sample types
fig7a <- ps %>%
  subset_samples(Pat_ID %in% focal_patients & day_cat!="Donor") %>%
  subset_taxa(Kingdom == "Bacteria") %>%
  aggregate_topn(level = "Species", detection = 0.01, topn = 15, toi = "Escherichia_coli") %>%
  myotudf(rankname = "Species", vars_of_interest = c("Pat_ID", "sample_type_spec", "day_cat")) %>%
  select(1:4, contains("Escherichia")) %>%
  mutate(Escherichia_coli = if_else(Escherichia_coli < 1e-6, 0, Escherichia_coli)) %>%
  ggplot(aes(x = day_cat, y = Escherichia_coli, color = Pat_ID, group = Pat_ID, shape = Pat_ID)) +
    geom_jitter(size = 3, width = 0.03, height = 0) +
    geom_line(linetype = 3, linewidth = 1) +
    scale_color_manual(values = patient_colors) +
    scale_y_continuous(
      name      = expression("Relative abundance " * italic("E. coli")),
      labels    = scales::label_percent(trim = FALSE),
      transform = "log10",
      limits    = c(1e-7, 1)
    ) +
    scale_x_discrete(name = NULL) +
    facet_grid(rows = "sample_type_spec", labeller = as_labeller(sample_type_labels)) +
    theme_minimal() +
    theme(
      axis.text       = element_text(size = 12),
      axis.title      = element_text(size = 12),
      strip.text      = element_text(size = 10),
      legend.text     = element_text(size = 10),
      legend.title    = element_blank(),
      legend.key.size = unit(5, "mm")
    )

# Fig. 9B: Bray-Curtis PCoA of recipient stool samples (Pat. 1, Pat. 2)
recipients <- c("Pat. 1", "Pat. 2")

pcoa <- ps %>%
  subset_samples(Pat_ID %in% recipients & sample_type_spec == "stool") %>%
  subset_taxa(Kingdom == "Bacteria") %>%
  mypcoa(norandom = TRUE, method = "PCoA", dist = "bray",
         test_group = "Pat_ID", shape = "day_cat", shape_keep_group = "day_cat",
         outputdir = outputdir)

fig7b <- pcoa[[1]]
fig7b$layers <- NULL
fig7b <- fig7b +
  geom_point(aes(shape = day_cat), size = 3, stroke = 1) +
  ggrepel::geom_text_repel(aes(label = day)) +
  scale_color_manual(values = patient_colors[recipients], breaks = recipients) +
  scale_shape_manual(values = c(19, 1, 6, 2, 5)) +
  labs(caption = NULL) +
  theme_minimal()

fig7a+fig7b
ggsave(filename = "./fig9.jpeg",width=297,height=210,units = "mm")

