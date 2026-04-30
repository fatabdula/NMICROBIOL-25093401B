# Extended Data 5 

pkgs <- c("readxl", "dplyr", "tidyr", "ggplot2", "stringr", "purrr", "scales", "tcltk")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(purrr)
library(scales)
library(tcltk)

file_path <- tclvalue(tkgetOpenFile(filetypes = "{{Excel Files} {.xlsx .xls}}"))
if (file_path == "") stop("No file selected")

save_dir <- tclvalue(tkchooseDirectory())
if (save_dir == "") stop("No save directory selected")

sheets <- excel_sheets(file_path)
print(data.frame(number = seq_along(sheets), sheet = sheets))
sheet_name <- sheets[as.integer(readline("Sheet number: "))]

clean_num <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "Missing", "missing", "nd", "ND")] <- NA
  suppressWarnings(as.numeric(x))
}

format_p <- function(x) {
  ifelse(is.na(x), NA_character_, formatC(x, format = "f", digits = 6))
}

safe_name <- function(x) {
  gsub("[^A-Za-z0-9_-]", "_", x)
}

raw <- read_excel(file_path, sheet = sheet_name)

needed <- c("Strain", "Serum", "PFU/mL")
missing <- setdiff(needed, names(raw))
if (length(missing) > 0) stop(paste("Missing columns:", paste(missing, collapse = ", ")))

serum_data <- raw %>%
  transmute(
    strain = as.character(Strain),
    serum = as.character(Serum),
    pfu_ml = clean_num(`PFU/mL`)
  ) %>%
  filter(!is.na(strain), !is.na(serum), !is.na(pfu_ml)) %>%
  mutate(
    patient = case_when(
      str_detect(serum, regex("Patient ?1|Patient1", ignore_case = TRUE)) ~ "Patient 1",
      str_detect(serum, regex("Patient ?2|Patient2", ignore_case = TRUE)) ~ "Patient 2",
      str_detect(serum, regex("Patient ?3|Patient3", ignore_case = TRUE)) ~ "Patient 3",
      TRUE ~ NA_character_
    ),
    visit = case_when(
      str_detect(serum, regex("vb|baseline|sot|start", ignore_case = TRUE)) ~ "vb",
      str_detect(serum, regex("v16|eot|end", ignore_case = TRUE)) ~ "v16",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(patient), !is.na(visit)) %>%
  mutate(
    patient = factor(patient, levels = c("Patient 1", "Patient 2", "Patient 3")),
    visit = factor(visit, levels = c("vb", "v16"))
  )

stats <- serum_data %>%
  group_by(strain, patient) %>%
  summarise(
    n_vb = sum(visit == "vb"),
    n_v16 = sum(visit == "v16"),
    mean_vb = mean(pfu_ml[visit == "vb"], na.rm = TRUE),
    mean_v16 = mean(pfu_ml[visit == "v16"], na.rm = TRUE),
    sd_vb = sd(pfu_ml[visit == "vb"], na.rm = TRUE),
    sd_v16 = sd(pfu_ml[visit == "v16"], na.rm = TRUE),
    p_value = tryCatch(
      t.test(
        pfu_ml[visit == "v16"],
        pfu_ml[visit == "vb"],
        alternative = "less",
        paired = FALSE,
        var.equal = FALSE
      )$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    test = "unpaired one-sided Welch t-test",
    comparison = "v16 < vb"
  )

plot_data <- serum_data %>%
  mutate(group = interaction(patient, visit, sep = " "))

make_plot <- function(selected_strain) {
  d <- plot_data %>% filter(strain == selected_strain)
  s <- stats %>% filter(strain == selected_strain)

  ggplot(d, aes(x = visit, y = pfu_ml)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.08, height = 0) +
    facet_wrap(~ patient, nrow = 1) +
    scale_y_continuous(trans = "log10", labels = scales::scientific) +
    labs(
      title = selected_strain,
      subtitle = paste("Welch t-test, one-sided: v16 < vb"),
      x = NULL,
      y = "PFU/mL"
    ) +
    theme_bw()
}

plots <- serum_data %>%
  distinct(strain) %>%
  pull(strain) %>%
  set_names() %>%
  map(make_plot)

for (nm in names(plots)) {
  print(plots[[nm]])
  ggsave(
    file.path(save_dir, paste0("ED5_serum_neutralization_", safe_name(nm), ".pdf")),
    plots[[nm]],
    width = 7,
    height = 4
  )
}

write.csv(serum_data, file.path(save_dir, "ED5_serum_neutralization_cleaned.csv"), row.names = FALSE)
write.csv(stats, file.path(save_dir, "ED5_serum_neutralization_statistics.csv"), row.names = FALSE)
print(stats %>% mutate(p_value = format_p(p_value)))

