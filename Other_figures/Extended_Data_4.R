# Extended Data 4 

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
names(raw) <- gsub("\\s+", " ", names(raw))

nm <- names(raw)

pick_col <- function(pattern) {
  hit <- grep(pattern, nm, ignore.case = TRUE, value = TRUE)
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

clean_num <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "Missing", "missing")] <- NA
  x[grepl("^\\d{1,3}(,\\d{3})+$", x)] <- gsub(",", "", x[grepl("^\\d{1,3}(,\\d{3})+$", x)])
  x[grepl("^\\d+,\\d+$", x)] <- sub(",", ".", x[grepl("^\\d+,\\d+$", x)])
  suppressWarnings(as.numeric(x))
}

patient_col <- nm[1]
day_col <- pick_col("^days$|^day$")
visit_col <- pick_col("^visit$")

parameter_cols <- c(
  pH = pick_col("^pH$"),
  leukocytes = pick_col("leuk"),
  erythrocytes = pick_col("ery"),
  bacteria = pick_col("bact")
)

parameter_cols <- parameter_cols[!is.na(parameter_cols)]

if (is.na(day_col)) stop("Day column was not found")
if (length(parameter_cols) == 0) stop("No urinalysis parameter columns were found")

urinalysis <- raw %>%
  mutate(row_id = row_number()) %>%
  transmute(
    patient = as.character(.data[[patient_col]]),
    day = clean_num(.data[[day_col]]),
    visit = if (!is.na(visit_col)) as.character(.data[[visit_col]]) else NA_character_,
    row_id,
    across(all_of(unname(parameter_cols)), clean_num)
  )

names(urinalysis)[match(unname(parameter_cols), names(urinalysis))] <- names(parameter_cols)

urinalysis_long <- urinalysis %>%
  mutate(
    patient = recode(
      patient,
      "Patient1" = "Patient 1",
      "Patient2" = "Patient 2",
      "Patient3" = "Patient 3",
      .default = patient
    )
  ) %>%
  pivot_longer(
    cols = all_of(names(parameter_cols)),
    names_to = "parameter",
    values_to = "value"
  ) %>%
  filter(!is.na(day), !is.na(value)) %>%
  group_by(patient, day, parameter) %>%
  arrange(row_id, .by_group = TRUE) %>%
  mutate(visit_order = row_number()) %>%
  ungroup() %>%
  mutate(
    value_plot = case_when(
      parameter == "bacteria" & value > 0 ~ log10(value),
      parameter == "bacteria" & value == 0 ~ 0,
      TRUE ~ value
    )
  )

timepoints <- urinalysis_long %>%
  distinct(day, visit_order) %>%
  arrange(day, visit_order) %>%
  mutate(
    timepoint = row_number(),
    label = paste0("V", visit_order)
  )

plot_data <- urinalysis_long %>%
  left_join(timepoints, by = c("day", "visit_order"))

make_line_plot <- function(df, selected_parameter) {
  d <- df %>% filter(parameter == selected_parameter)

  y_label <- if (selected_parameter == "bacteria") {
    "Bacteria, log10 scale; zeros kept as 0"
  } else {
    selected_parameter
  }

  ggplot(d, aes(timepoint, value_plot, group = patient, colour = patient)) +
    geom_line() +
    geom_point() +
    scale_x_continuous(
      breaks = timepoints$timepoint,
      labels = paste0("Day ", timepoints$day, " ", timepoints$label)
    ) +
    labs(x = NULL, y = y_label, colour = "Patient") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plots <- lapply(unique(plot_data$parameter), function(x) make_line_plot(plot_data, x))
names(plots) <- unique(plot_data$parameter)

for (nm in names(plots)) {
  print(plots[[nm]])
  ggsave(
    file.path(save_dir, paste0("ED4_", nm, "_line_plot.pdf")),
    plots[[nm]],
    width = 7,
    height = 4
  )
}

write.csv(plot_data, file.path(save_dir, "ED4_urinalysis_long.csv"), row.names = FALSE)

