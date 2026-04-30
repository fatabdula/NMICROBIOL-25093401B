# Extended Data 3 

pkgs <- c("readxl", "dplyr", "tidyr", "ggplot2", "stringr", "tcltk")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

library(readxl)
library(dplyr)
library(tidyr)
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

names(raw) <- clean_names(names(raw))
nm <- names(raw)

pick_col <- function(keys) {
  hit <- unlist(lapply(keys, function(k) nm[grepl(k, nm, ignore.case = TRUE)]))
  hit <- unique(hit)
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

patient_col <- nm[1]
day_col <- pick_col(c("^days$", "^day$"))
visit_col <- pick_col(c("^visit$"))
time_col <- pick_col(c("^time$"))
symptom_col <- pick_col(c("^s_total$", "symptomatic"))
dd_col <- pick_col(c("^dd_total$", "differential"))
qol_col <- pick_col(c("^qol_total$", "quality"))

needed <- c(day_col, visit_col, time_col, symptom_col, dd_col, qol_col)
if (any(is.na(needed))) stop("Required columns were not found")

clock_minutes <- function(x) {
  x <- trimws(as.character(x))
  out <- rep(NA_real_, length(x))

  ok <- grepl("^\\d{1,2}:\\d{2}(:\\d{2})?$", x)
  z <- str_split_fixed(x[ok], ":", 3)

  h <- suppressWarnings(as.numeric(z[, 1]))
  m <- suppressWarnings(as.numeric(z[, 2]))
  s <- suppressWarnings(as.numeric(z[, 3]))
  s[is.na(s)] <- 0

  out[ok] <- h * 60 + m + s / 60
  out
}

acss_long <- raw %>%
  transmute(
    patient = as.character(.data[[patient_col]]),
    day = as.numeric(.data[[day_col]]),
    visit = as.numeric(.data[[visit_col]]),
    time = as.character(.data[[time_col]]),
    symptomatic = as.numeric(.data[[symptom_col]]),
    differential_diagnosis = as.numeric(.data[[dd_col]]),
    quality_of_life = as.numeric(.data[[qol_col]])
  ) %>%
  mutate(
    patient = recode(patient, "Patient1" = "Patient 1", "Patient2" = "Patient 2", "Patient3" = "Patient 3"),
    time_min = clock_minutes(time)
  ) %>%
  arrange(day, visit, time_min, patient) %>%
  mutate(timepoint = paste0("Day ", day, " / Visit ", visit)) %>%
  pivot_longer(
    c(symptomatic, differential_diagnosis, quality_of_life),
    names_to = "score_type",
    values_to = "score"
  ) %>%
  mutate(
    score = replace_na(score, 0),
    score_type = factor(
      score_type,
      levels = c("symptomatic", "differential_diagnosis", "quality_of_life")
    )
  )

timepoint_order <- acss_long %>%
  distinct(day, visit, time_min, timepoint) %>%
  arrange(day, visit, time_min) %>%
  pull(timepoint)

acss_long <- acss_long %>%
  mutate(
    timepoint = factor(timepoint, levels = unique(timepoint_order)),
    patient = factor(patient, levels = c("Patient 1", "Patient 2", "Patient 3"))
  )

p_acss <- ggplot(acss_long, aes(timepoint, score, fill = score_type)) +
  geom_col() +
  facet_wrap(~ patient, ncol = 1) +
  labs(x = NULL, y = "ACSS score", fill = "Score type") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_acss)

write.csv(acss_long, file.path(save_dir, "ED3_ACSS_long.csv"), row.names = FALSE)
ggsave(file.path(save_dir, "ED3_ACSS_bar_plot.pdf"), p_acss, width = 8, height = 6)

