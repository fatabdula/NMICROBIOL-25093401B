# Figure 2 

pkgs <- c("readxl", "dplyr", "tidyr", "ggplot2", "tcltk")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tcltk)

file_path <- tclvalue(tkgetOpenFile(filetypes = "{{Excel Files} {.xlsx .xls}}"))
if (file_path == "") stop("No file selected")

save_dir <- tclvalue(tkchooseDirectory())
if (save_dir == "") stop("No save directory selected")

sheets <- excel_sheets(file_path)
print(data.frame(number = seq_along(sheets), sheet = sheets))

sheet_ab <- sheets[as.integer(readline("Sheet for Figure 2a-b: "))]
sheet_c  <- sheets[as.integer(readline("Sheet for Figure 2c: "))]

clean_names <- function(x) {
  tolower(gsub("[^a-z0-9]+", "_", x))
}

pick_col <- function(nm, keys) {
  hit <- unlist(lapply(keys, function(k) nm[grepl(k, nm, ignore.case = TRUE)]))
  hit <- unique(hit)
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

visit_offset <- function(n) {
  if (n == 1) return(0)
  seq(-0.18, 0.18, length.out = n)
}

raw_ab <- read_excel(file_path, sheet = sheet_ab)
names(raw_ab) <- clean_names(names(raw_ab))

nm <- names(raw_ab)

sample_col  <- pick_col(nm, c("^sample_name$"))
case_col    <- pick_col(nm, c("^case$"))
patient_col <- pick_col(nm, c("^patient$"))
days_col    <- pick_col(nm, c("^days$", "^day$"))
visit_col   <- pick_col(nm, c("^visit$"))
phage_col   <- pick_col(nm, c("results_phage", "phage_total", "phage"))
bact_col    <- pick_col(nm, c("results_bacteria", "bacteria_cfu", "cfu"))

needed_ab <- c(days_col, phage_col, bact_col)
if (any(is.na(needed_ab))) stop("Required Figure 2a-b columns were not found")

fig2ab <- raw_ab %>%
  mutate(
    patient = coalesce(
      if (!is.na(patient_col)) as.character(.data[[patient_col]]) else NA_character_,
      if (!is.na(case_col)) as.character(.data[[case_col]]) else NA_character_,
      if (!is.na(sample_col)) as.character(.data[[sample_col]]) else NA_character_
    ),
    day = as.numeric(.data[[days_col]]),
    visit = if (!is.na(visit_col)) as.numeric(.data[[visit_col]]) else NA_real_,
    phage_total = as.numeric(.data[[phage_col]]),
    ecoli_cfu_ml = as.numeric(gsub(",", "", as.character(.data[[bact_col]])))
  ) %>%
  mutate(
    patient = gsub(".*(Patient ?[123]).*", "\\1", patient),
    patient = gsub("Patient([123])", "Patient \\1", patient),
    patient = recode(
      patient,
      "Patient1" = "Patient 1",
      "Patient2" = "Patient 2",
      "Patient3" = "Patient 3",
      .default = patient
    )
  ) %>%
  filter(!is.na(patient), !is.na(day)) %>%
  group_by(patient, day) %>%
  arrange(visit, .by_group = TRUE) %>%
  mutate(
    visit_order = row_number(),
    x_day = day + visit_offset(n())[visit_order]
  ) %>%
  ungroup() %>%
  mutate(
    ecoli_log10 = ifelse(ecoli_cfu_ml > 0, log10(ecoli_cfu_ml), 0),
    phage_log10 = ifelse(phage_total > 0, log10(phage_total), 0)
  )

p_ecoli <- ggplot(fig2ab, aes(x_day, ecoli_log10, group = patient, colour = patient)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = sort(unique(fig2ab$day))) +
  labs(x = "Days", y = "E. coli, log10 CFU/mL") +
  theme_bw()

p_phage <- ggplot(fig2ab, aes(x_day, phage_log10, group = patient, colour = patient)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = sort(unique(fig2ab$day))) +
  labs(x = "Days", y = "Phage, log10 total PFU") +
  theme_bw()

raw_c <- read_excel(file_path, sheet = sheet_c)
names(raw_c) <- clean_names(names(raw_c))

make_uti_long <- function(df) {
  nm <- names(df)

  if (all(c("patient", "pre", "post") %in% nm)) {
    df %>%
      select(patient, pre, post) %>%
      pivot_longer(c(pre, post), names_to = "period", values_to = "uti_per_year")
  } else {
    y_col <- intersect(c("utis_per_year", "utis_year", "utis"), nm)[1]

    if (!all(c("patient", "treatment") %in% nm) || is.na(y_col)) {
      stop("Required Figure 2c columns were not found")
    }

    df %>%
      transmute(
        patient = .data[["patient"]],
        period = .data[["treatment"]],
        uti_per_year = .data[[y_col]]
      )
  }
}

fig2c <- make_uti_long(raw_c) %>%
  mutate(
    patient = as.character(patient),
    patient = gsub(".*(Patient ?[123]).*", "\\1", patient),
    patient = gsub("Patient([123])", "Patient \\1", patient),
    patient = recode(
      patient,
      "1" = "Patient 1",
      "2" = "Patient 2",
      "3" = "Patient 3",
      "Patient1" = "Patient 1",
      "Patient2" = "Patient 2",
      "Patient3" = "Patient 3",
      .default = patient
    ),
    period = factor(tolower(as.character(period)), levels = c("pre", "post"), labels = c("Pre", "Post")),
    uti_per_year = as.numeric(uti_per_year)
  ) %>%
  filter(!is.na(patient), !is.na(period), !is.na(uti_per_year))

p_uti <- ggplot(fig2c, aes(period, uti_per_year, group = patient, colour = patient)) +
  geom_line() +
  geom_point(size = 3) +
  labs(x = "Treatment period", y = "UTIs per year") +
  theme_bw()

print(p_ecoli)
print(p_phage)
print(p_uti)

write.csv(fig2ab, file.path(save_dir, "Figure_2ab_processed.csv"), row.names = FALSE)
write.csv(fig2c, file.path(save_dir, "Figure_2c_processed.csv"), row.names = FALSE)

ggsave(file.path(save_dir, "Figure_2a_ecoli_line_plot.pdf"), p_ecoli, width = 6, height = 4)
ggsave(file.path(save_dir, "Figure_2b_phage_line_plot.pdf"), p_phage, width = 6, height = 4)
ggsave(file.path(save_dir, "Figure_2c_uti_line_plot.pdf"), p_uti, width = 5, height = 4)
