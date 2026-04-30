# Extended Data 6 

pkgs <- c("readxl", "dplyr", "ggplot2", "stringr", "tcltk")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

library(readxl)
library(dplyr)
library(ggplot2)
library(stringr)
library(tcltk)

file_path <- tclvalue(tkgetOpenFile(filetypes = "{{Excel Files} {.xlsx .xls}}"))
if (file_path == "") stop("No file selected")

save_dir <- tclvalue(tkchooseDirectory())
if (save_dir == "") stop("No save directory selected")

sheets <- excel_sheets(file_path)
print(data.frame(number = seq_along(sheets), sheet = sheets))
sheet_name <- sheets[as.integer(readline("Sheet number: "))]

raw <- read_excel(file_path, sheet = sheet_name)

clean_names <- function(x) {
  tolower(gsub("[^a-z0-9]+", "_", x))
}

clean_num <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "Missing", "missing")] <- NA
  suppressWarnings(as.numeric(x))
}

names(raw) <- clean_names(names(raw))
nm <- names(raw)

pick_col <- function(keys) {
  hit <- unlist(lapply(keys, function(k) nm[grepl(k, nm, ignore.case = TRUE)]))
  hit <- unique(hit)
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

sample_col <- pick_col(c("sample_name", "sample"))
hour_col <- pick_col(c("^hour$", "time"))
phage_col <- pick_col(c("results_phage", "phage_total", "phage"))

needed <- c(sample_col, hour_col, phage_col)
if (any(is.na(needed))) stop("Required columns were not found")

phage_recovery <- raw %>%
  transmute(
    sample = as.character(.data[[sample_col]]),
    hour = clean_num(.data[[hour_col]]),
    phage_total = clean_num(.data[[phage_col]])
  ) %>%
  filter(!is.na(sample), !is.na(hour), !is.na(phage_total)) %>%
  mutate(
    patient = case_when(
      str_detect(sample, regex("Patient ?1|Patient1", ignore_case = TRUE)) ~ "Patient 1",
      str_detect(sample, regex("Patient ?2|Patient2", ignore_case = TRUE)) ~ "Patient 2",
      str_detect(sample, regex("Patient ?3|Patient3", ignore_case = TRUE)) ~ "Patient 3",
      TRUE ~ NA_character_
    ),
    patient = factor(patient, levels = c("Patient 1", "Patient 2", "Patient 3")),
    hour = factor(hour, levels = sort(unique(hour))),
    phage_log10 = case_when(
      phage_total > 0 ~ log10(phage_total),
      phage_total == 0 ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(patient), !is.na(phage_log10))

p_recovery <- ggplot(phage_recovery, aes(hour, phage_log10, fill = patient)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  labs(x = "Hours", y = "Phage, log10 total PFU", fill = "Patient") +
  theme_bw()

print(p_recovery)

write.csv(phage_recovery, file.path(save_dir, "ED6_phage_recovery.csv"), row.names = FALSE)
ggsave(file.path(save_dir, "ED6_phage_recovery_bar_plot.pdf"), p_recovery, width = 6, height = 4)

