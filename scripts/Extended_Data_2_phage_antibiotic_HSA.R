# Extended Data 2 

pkgs <- c("readxl", "dplyr", "tidyr", "stringr", "ggplot2", "tcltk")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(tcltk)

raw_file <- tclvalue(tkgetOpenFile(
  title = "Select raw OD600 workbook",
  filetypes = "{{Excel Files} {.xlsx .xls}}"
))
if (raw_file == "") stop("No raw OD600 workbook selected")

layout_file <- tclvalue(tkgetOpenFile(
  title = "Select plate-layout workbook",
  filetypes = "{{Excel Files} {.xlsx .xls}}"
))
if (layout_file == "") stop("No plate-layout workbook selected")

save_dir <- tclvalue(tkchooseDirectory(title = "Select output folder"))
if (save_dir == "") stop("No output folder selected")

raw_sheets <- excel_sheets(raw_file)
layout_sheets <- excel_sheets(layout_file)

print(data.frame(number = seq_along(raw_sheets), raw_sheet = raw_sheets))
raw_sheet <- raw_sheets[as.integer(readline("Raw OD600 sheet number: "))]

print(data.frame(number = seq_along(layout_sheets), layout_sheet = layout_sheets))
layout_sheet <- layout_sheets[as.integer(readline("Layout sheet number: "))]

auc_window_h <- as.numeric(readline("AUC window in hours, e.g. 24: "))
if (is.na(auc_window_h)) auc_window_h <- 24

parse_time_h <- function(x) {
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

auc_trapz <- function(time_h, y, end_h) {
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

first_number <- function(x) {
  x <- str_replace_all(as.character(x), ",", ".")
  suppressWarnings(as.numeric(str_extract(x, "[0-9]+(?:\\.[0-9]+)?(?:[eE][+-]?[0-9]+)?")))
}

read_od_sheet <- function(path, sheet) {
  preview <- suppressMessages(read_excel(path, sheet = sheet, n_max = 50, col_names = FALSE))

  header_row <- which(apply(preview, 1, function(z) {
    txt <- paste(z, collapse = " ")
    str_detect(txt, regex("Time", ignore_case = TRUE)) && str_detect(txt, regex("A1", ignore_case = TRUE))
  }))[1]

  if (is.na(header_row)) header_row <- 1

  dat <- suppressMessages(read_excel(path, sheet = sheet, skip = header_row - 1, .name_repair = "unique"))
  names(dat) <- make.unique(names(dat))

  time_col <- names(dat)[str_detect(names(dat), regex("^time$", ignore_case = TRUE))][1]
  names(dat)[names(dat) == time_col] <- "Time"

  well_cols <- names(dat)[str_detect(names(dat), "^[A-Ha-h][0-9]+$")]
  names(dat)[match(well_cols, names(dat))] <- toupper(well_cols)

  dat %>% select(Time, matches("^[A-H][0-9]+$"))
}

read_layout <- function(path, sheet) {
  lay <- suppressMessages(read_excel(path, sheet = sheet, col_names = FALSE))
  lay <- as.data.frame(lay, stringsAsFactors = FALSE)

  row_id <- which(as.character(lay[[1]]) %in% LETTERS[1:8])
  if (length(row_id) != 8) stop("Could not find rows A-H in the layout sheet")

  out <- list()

  for (i in row_id) {
    row_name <- as.character(lay[i, 1])

    for (j in 1:12) {
      out[[length(out) + 1]] <- data.frame(
        well = paste0(row_name, j),
        row = row_name,
        column = j,
        layout_text = as.character(lay[i, j + 1]),
        stringsAsFactors = FALSE
      )
    }
  }

  bind_rows(out)
}

layout <- read_layout(layout_file, layout_sheet) %>%
  mutate(
    condition = case_when(
      str_detect(layout_text, regex("blank|media|LB", ignore_case = TRUE)) ~ "blank",
      str_detect(layout_text, regex("GC|growth", ignore_case = TRUE)) ~ "growth_control",
      row %in% LETTERS[1:7] & column <= 10 ~ "combo",
      row %in% LETTERS[1:7] & column == 11 ~ "phage_only",
      row == "H" & column <= 10 ~ "antibiotic_only",
      TRUE ~ "other"
    ),
    antibiotic_ug_ml = if_else(condition %in% c("combo", "antibiotic_only"), first_number(layout_text), NA_real_),
    phage_pfu_ml = if_else(condition %in% c("combo", "phage_only"), first_number(layout_text), NA_real_)
  )

od_long <- read_od_sheet(raw_file, raw_sheet) %>%
  mutate(time_h = parse_time_h(Time)) %>%
  select(time_h, matches("^[A-H][0-9]+$")) %>%
  pivot_longer(-time_h, names_to = "well", values_to = "od600") %>%
  mutate(od600 = as.numeric(od600)) %>%
  left_join(layout, by = "well") %>%
  filter(!is.na(time_h), !is.na(od600), !is.na(condition))

growth_curves <- od_long %>%
  group_by(condition, antibiotic_ug_ml, phage_pfu_ml, time_h) %>%
  summarise(mean_od600 = mean(od600, na.rm = TRUE), .groups = "drop")

auc_table <- od_long %>%
  group_by(well, row, column, condition, antibiotic_ug_ml, phage_pfu_ml, layout_text) %>%
  summarise(auc = auc_trapz(time_h, od600, auc_window_h), .groups = "drop")

blank_ref <- median(auc_table$auc[auc_table$condition == "blank"], na.rm = TRUE)
growth_ref <- median(auc_table$auc[auc_table$condition == "growth_control"], na.rm = TRUE)

auc_table <- auc_table %>%
  mutate(
    pi = 100 * (1 - (auc - blank_ref) / (growth_ref - blank_ref)),
    pi = pmin(100, pmax(0, pi))
  )

phage_only <- auc_table %>%
  filter(condition == "phage_only") %>%
  select(row, auc_phage = auc, pi_phage = pi)

antibiotic_only <- auc_table %>%
  filter(condition == "antibiotic_only") %>%
  select(column, auc_antibiotic = auc, pi_antibiotic = pi)

hsa_table <- auc_table %>%
  filter(condition == "combo") %>%
  left_join(phage_only, by = "row") %>%
  left_join(antibiotic_only, by = "column") %>%
  mutate(
    best_single_auc = pmin(auc_phage, auc_antibiotic, na.rm = TRUE),
    best_single_auc = ifelse(is.finite(best_single_auc), best_single_auc, NA_real_),
    delta_auc = best_single_auc - auc,
    hsa_threshold = abs(best_single_auc) * 0.10,
    hsa_class = case_when(
      !is.finite(delta_auc) ~ NA_character_,
      delta_auc > hsa_threshold ~ "Additive",
      delta_auc < -hsa_threshold ~ "Antagonistic",
      TRUE ~ "No difference"
    )
  )

hsa_summary <- hsa_table %>%
  filter(!is.na(hsa_class)) %>%
  count(hsa_class) %>%
  mutate(percent = 100 * n / sum(n))

p_growth <- ggplot(growth_curves, aes(time_h, mean_od600, group = interaction(condition, antibiotic_ug_ml, phage_pfu_ml), colour = condition)) +
  geom_line() +
  labs(x = "Time [h]", y = "OD600") +
  theme_bw()

p_hsa <- ggplot(hsa_summary, aes(hsa_class, percent)) +
  geom_col() +
  geom_text(aes(label = paste0(n, " (", round(percent, 1), "%)")), vjust = -0.3) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = NULL, y = "Percent") +
  theme_bw()

print(p_growth)
print(p_hsa)

write.csv(od_long, file.path(save_dir, "ED2_OD600_long.csv"), row.names = FALSE)
write.csv(growth_curves, file.path(save_dir, "ED2_growth_curves.csv"), row.names = FALSE)
write.csv(auc_table, file.path(save_dir, "ED2_AUC_PI_table.csv"), row.names = FALSE)
write.csv(hsa_table, file.path(save_dir, "ED2_HSA_table.csv"), row.names = FALSE)
write.csv(hsa_summary, file.path(save_dir, "ED2_HSA_summary.csv"), row.names = FALSE)

ggsave(file.path(save_dir, "ED2_growth_curves.pdf"), p_growth, width = 8, height = 5)
ggsave(file.path(save_dir, "ED2_HSA_summary.pdf"), p_hsa, width = 5, height = 4)

