# Extended Data 1 

pkgs <- c("readxl", "dplyr", "tidyr", "ggplot2", "hms", "tcltk")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(hms)
library(tcltk)

file_path <- tclvalue(tkgetOpenFile(filetypes = "{{Excel Files} {.xlsx .xls}}"))
if (file_path == "") stop("No file selected")

save_dir <- tclvalue(tkchooseDirectory())
if (save_dir == "") stop("No save directory selected")

sheets <- excel_sheets(file_path)
print(data.frame(number = seq_along(sheets), sheet = sheets))
sheet_name <- sheets[as.integer(readline("Sheet number: "))]

raw <- read_excel(file_path, sheet = sheet_name)
well_cols <- names(raw)[grepl("^[A-F](?:[1-9]|1[0-2])$", names(raw))]

parse_time_hours <- function(x) {
  if (inherits(x, "hms")) return(as.numeric(x) / 3600)
  if (is.numeric(x)) return(as.numeric(x) * 24)

  x <- trimws(as.character(x))
  out <- suppressWarnings(as.numeric(x) * 24)
  hms_id <- grepl("^\\d{1,2}:\\d{2}:\\d{2}$", x)

  if (any(hms_id)) {
    out[hms_id] <- as.numeric(hms::as_hms(x[hms_id])) / 3600
  }

  out
}

auc_trap <- function(time, y, end_time) {
  keep <- !is.na(time) & !is.na(y)
  time <- time[keep]
  y <- y[keep]

  if (length(time) < 2) return(NA_real_)

  o <- order(time)
  time <- time[o]
  y <- y[o]

  if (max(time) < end_time) return(NA_real_)

  y_end <- approx(time, y, xout = end_time, ties = "ordered")$y
  keep <- time <= end_time

  time <- c(time[keep], end_time)
  y <- c(y[keep], y_end)

  id <- !duplicated(time)
  time <- time[id]
  y <- y[id]

  sum(diff(time) * (head(y, -1) + tail(y, -1)) / 2)
}

tra_long <- raw %>%
  mutate(time_h = parse_time_hours(Time)) %>%
  mutate(across(all_of(well_cols), ~ as.numeric(as.character(.x)))) %>%
  filter(!is.na(time_h), time_h <= 18) %>%
  select(time_h, all_of(well_cols)) %>%
  pivot_longer(all_of(well_cols), names_to = "well", values_to = "od600") %>%
  mutate(
    row_block = sub("([A-F]).*", "\\1", well),
    column = as.integer(sub("[A-F]", "", well)),
    treatment = case_when(
      column %in% 1:3 ~ "WT cocktail",
      column %in% 4:6 ~ "uroCOLE7-01",
      column %in% 7:9 ~ "Growth control",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(od600), !is.na(treatment))

growth_summary <- tra_long %>%
  group_by(time_h, row_block, treatment) %>%
  summarise(
    mean_od600 = mean(od600, na.rm = TRUE),
    sd_od600 = sd(od600, na.rm = TRUE),
    .groups = "drop"
  )

auc_table <- growth_summary %>%
  group_by(row_block, treatment) %>%
  summarise(
    auc_5h = auc_trap(time_h, mean_od600, 5),
    auc_18h = auc_trap(time_h, mean_od600, 18),
    .groups = "drop"
  )

growth_control <- auc_table %>%
  filter(treatment == "Growth control") %>%
  select(row_block, auc_5h_gc = auc_5h, auc_18h_gc = auc_18h)

pi_table <- auc_table %>%
  filter(treatment != "Growth control") %>%
  left_join(growth_control, by = "row_block") %>%
  mutate(
    pi_5h = 100 * (auc_5h_gc - auc_5h) / auc_5h_gc,
    pi_18h = 100 * (auc_18h_gc - auc_18h) / auc_18h_gc,
    pi_5h = pmax(pi_5h, 0),
    pi_18h = pmax(pi_18h, 0)
  ) %>%
  pivot_longer(c(pi_5h, pi_18h), names_to = "timepoint", values_to = "pi") %>%
  mutate(timepoint = recode(timepoint, pi_5h = "5 h", pi_18h = "18 h"))

p_growth <- ggplot(growth_summary, aes(time_h, mean_od600, colour = treatment, group = treatment)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ row_block) +
  labs(x = "Time [h]", y = "OD600") +
  theme_bw()

p_pi <- ggplot(pi_table, aes(row_block, pi, fill = treatment)) +
  geom_col(position = "dodge") +
  facet_wrap(~ timepoint) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = "Block", y = "PI [%]") +
  theme_bw()

print(p_growth)
print(p_pi)

write.csv(growth_summary, file.path(save_dir, "ED1_growth_summary.csv"), row.names = FALSE)
write.csv(auc_table, file.path(save_dir, "ED1_auc_table.csv"), row.names = FALSE)
write.csv(pi_table, file.path(save_dir, "ED1_pi_table.csv"), row.names = FALSE)

ggsave(file.path(save_dir, "ED1_growth_curves.pdf"), p_growth, width = 9, height = 5)
ggsave(file.path(save_dir, "ED1_PI_bar_plot.pdf"), p_pi, width = 7, height = 4)

