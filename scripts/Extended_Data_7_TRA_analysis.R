# Extended Data 7 

pkgs <- c("readxl", "dplyr", "tidyr", "ggplot2", "stringr", "hms", "tcltk")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(hms)
library(tcltk)

file_path <- tclvalue(tkgetOpenFile(filetypes = "{{Excel Files} {.xlsx .xls}}"))
if (file_path == "") stop("No file selected")

save_dir <- tclvalue(tkchooseDirectory())
if (save_dir == "") stop("No save directory selected")

sheets <- excel_sheets(file_path)
print(data.frame(number = seq_along(sheets), sheet = sheets))
selected <- readline("Sheet numbers, separated by comma: ")
sheet_id <- as.integer(strsplit(selected, ",")[[1]])
sheet_names <- sheets[sheet_id]

auc_windows <- c(5, 18)

treatment_layout <- function() {
  expand.grid(row = LETTERS[1:8], column = 1:12) %>%
    mutate(
      well = paste0(row, column),
      strain = paste0("Strain ", row),
      treatment = case_when(
        column %in% 1:3 ~ "Control",
        column %in% 4:6 ~ "E2",
        column %in% 7:9 ~ "phi41S",
        column %in% 10:12 ~ "Mix",
        TRUE ~ NA_character_
      ),
      replicate = ((column - 1) %% 3) + 1
    )
}

parse_time_h <- function(x) {
  if (inherits(x, "hms")) return(as.numeric(x) / 3600)
  if (inherits(x, "difftime")) return(as.numeric(x, units = "hours"))

  if (is.numeric(x)) {
    if (max(x, na.rm = TRUE) <= 2) return(x * 24)
    return(x)
  }

  x <- as.character(x)

  if (any(str_detect(x, ":"), na.rm = TRUE)) {
    z <- str_split(x, ":", simplify = TRUE)
    h <- suppressWarnings(as.numeric(z[, 1]))
    m <- suppressWarnings(as.numeric(z[, 2]))
    s <- suppressWarnings(as.numeric(z[, 3]))
    m[is.na(m)] <- 0
    s[is.na(s)] <- 0
    return(h + m / 60 + s / 3600)
  }

  suppressWarnings(as.numeric(x))
}

auc_trap <- function(time_h, y, end_h) {
  keep <- is.finite(time_h) & is.finite(y)
  time_h <- time_h[keep]
  y <- y[keep]

  if (length(time_h) < 2) return(NA_real_)

  o <- order(time_h)
  time_h <- time_h[o]
  y <- y[o]
  time_h <- time_h - min(time_h)

  if (max(time_h) < end_h) return(NA_real_)

  y_end <- approx(time_h, y, xout = end_h, ties = "ordered")$y
  keep <- time_h <= end_h

  time_h <- c(time_h[keep], end_h)
  y <- c(y[keep], y_end)

  id <- !duplicated(time_h)
  time_h <- time_h[id]
  y <- y[id]

  sum(diff(time_h) * (head(y, -1) + tail(y, -1)) / 2)
}

read_od_sheet <- function(path, sheet) {
  raw <- read_excel(path, sheet = sheet)
  names(raw) <- make.unique(names(raw))

  time_col <- names(raw)[str_detect(names(raw), regex("^time$", ignore_case = TRUE))][1]
  if (is.na(time_col)) stop(paste("No Time column found in", sheet))

  names(raw)[names(raw) == time_col] <- "Time"

  well_cols <- names(raw)[str_detect(names(raw), "^[A-Ha-h][0-9]+$")]
  names(raw)[match(well_cols, names(raw))] <- toupper(well_cols)

  raw %>% select(Time, matches("^[A-H][0-9]+$"))
}

layout <- treatment_layout()

od_long <- lapply(sheet_names, function(sh) {
  read_od_sheet(file_path, sh) %>%
    mutate(time_h = parse_time_h(Time)) %>%
    select(time_h, matches("^[A-H][0-9]+$")) %>%
    pivot_longer(-time_h, names_to = "well", values_to = "od600") %>%
    mutate(
      plate = sh,
      od600 = as.numeric(od600)
    ) %>%
    left_join(layout, by = "well") %>%
    filter(!is.na(time_h), !is.na(od600), !is.na(treatment))
}) %>%
  bind_rows()

growth_curves <- od_long %>%
  group_by(plate, strain, treatment, time_h) %>%
  summarise(
    mean_od600 = mean(od600, na.rm = TRUE),
    sd_od600 = sd(od600, na.rm = TRUE),
    .groups = "drop"
  )

auc_table <- lapply(auc_windows, function(h) {
  od_long %>%
    group_by(plate, strain, row, treatment, replicate, well) %>%
    summarise(auc = auc_trap(time_h, od600, h), .groups = "drop") %>%
    mutate(auc_window_h = h)
}) %>%
  bind_rows()

control_auc <- auc_table %>%
  filter(treatment == "Control") %>%
  group_by(plate, strain, row, auc_window_h) %>%
  summarise(control_auc = median(auc, na.rm = TRUE), .groups = "drop")

pi_table <- auc_table %>%
  left_join(control_auc, by = c("plate", "strain", "row", "auc_window_h")) %>%
  mutate(
    pi = ifelse(control_auc > 0, 100 * (1 - auc / control_auc), NA_real_),
    pi = pmin(100, pmax(0, pi))
  )

pi_summary <- pi_table %>%
  filter(treatment != "Control") %>%
  group_by(plate, strain, treatment, auc_window_h) %>%
  summarise(
    pi_mean = mean(pi, na.rm = TRUE),
    pi_sd = sd(pi, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

p_growth <- ggplot(growth_curves, aes(time_h, mean_od600, colour = treatment, group = treatment)) +
  geom_line() +
  geom_point() +
  facet_grid(plate ~ strain) +
  labs(x = "Time [h]", y = "OD600") +
  theme_bw()

p_pi <- ggplot(pi_summary, aes(treatment, pi_mean, fill = treatment)) +
  geom_col() +
  geom_errorbar(aes(ymin = pmax(0, pi_mean - pi_sd), ymax = pmin(100, pi_mean + pi_sd)), width = 0.2) +
  facet_grid(auc_window_h + plate ~ strain, labeller = label_both) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = NULL, y = "PI [%]") +
  theme_bw()

print(p_growth)
print(p_pi)

write.csv(od_long, file.path(save_dir, "ED7_OD600_long.csv"), row.names = FALSE)
write.csv(growth_curves, file.path(save_dir, "ED7_growth_curves.csv"), row.names = FALSE)
write.csv(auc_table, file.path(save_dir, "ED7_AUC_table.csv"), row.names = FALSE)
write.csv(pi_table, file.path(save_dir, "ED7_PI_table.csv"), row.names = FALSE)
write.csv(pi_summary, file.path(save_dir, "ED7_PI_summary.csv"), row.names = FALSE)

ggsave(file.path(save_dir, "ED7_growth_curves.pdf"), p_growth, width = 10, height = 5)
ggsave(file.path(save_dir, "ED7_PI_bar_plot.pdf"), p_pi, width = 10, height = 5)

