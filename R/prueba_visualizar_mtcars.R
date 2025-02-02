library("sessioninfo")
library("here")
library("ggplot2")

## Hello world
print("Soy Silvana")

## Directorios
dir_plots <- here::here("figuras")
dir_rdata <- here::here("processed-data")

## Crear directorio para las figuras y archivos
dir.create(dir_plots, showWarnings = FALSE)
dir.create(dir_rdata, showWarnings = FALSE)

## Hacer una imagen de ejemplo
pdf(file.path(dir_plots, "mtcars_gear_vs_mpg.pdf"),
    useDingbats = FALSE
)
ggplot(mtcars, aes(group = gear, y = mpg)) +
  geom_boxplot()
dev.off()

## Para reproducir mi código
options(width = 120)
sessioninfo::session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.4.2 (2024-10-31)
# os       Ubuntu 20.04.3 LTS
# system   x86_64, linux-gnu
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       America/Mexico_City
# date     2025-02-01
# rstudio  2024.12.0+467 Kousa Dogwood (server)
# pandoc   3.2 @ /usr/lib/rstudio-server/bin/quarto/bin/tools/x86_64/ (via rmarkdown)
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package     * version date (UTC) lib source
# callr         3.7.6   2024-03-25 [2] CRAN (R 4.4.2)
# cli           3.6.3   2024-06-21 [2] CRAN (R 4.4.2)
# colorspace    2.1-1   2024-07-26 [2] CRAN (R 4.4.2)
# digest        0.6.37  2024-08-19 [2] CRAN (R 4.4.2)
# dplyr         1.1.4   2023-11-17 [2] CRAN (R 4.4.2)
# evaluate      1.0.3   2025-01-10 [2] CRAN (R 4.4.2)
# farver        2.1.2   2024-05-13 [2] CRAN (R 4.4.2)
# fastmap       1.2.0   2024-05-15 [2] CRAN (R 4.4.2)
# fs            1.6.5   2024-10-30 [2] CRAN (R 4.4.2)
# generics      0.1.3   2022-07-05 [2] CRAN (R 4.1.2)
# ggplot2     * 3.5.1   2024-04-23 [2] CRAN (R 4.4.2)
# glue          1.8.0   2024-09-30 [2] CRAN (R 4.4.2)
# gtable        0.3.6   2024-10-25 [2] CRAN (R 4.4.2)
# here        * 1.0.1   2020-12-13 [2] CRAN (R 4.1.2)
# htmltools     0.5.8.1 2024-04-04 [2] CRAN (R 4.4.2)
# knitr         1.49    2024-11-08 [2] CRAN (R 4.4.2)
# labeling      0.4.3   2023-08-29 [2] CRAN (R 4.4.2)
# lifecycle     1.0.4   2023-11-07 [2] CRAN (R 4.4.2)
# magrittr      2.0.3   2022-03-30 [2] CRAN (R 4.4.2)
# munsell       0.5.1   2024-04-01 [2] CRAN (R 4.4.2)
# pillar        1.10.1  2025-01-07 [2] CRAN (R 4.4.2)
# pkgconfig     2.0.3   2019-09-22 [2] CRAN (R 4.0.3)
# processx      3.8.5   2025-01-08 [2] CRAN (R 4.4.2)
# ps            1.8.1   2024-10-28 [2] CRAN (R 4.4.2)
# purrr         1.0.2   2023-08-10 [2] CRAN (R 4.4.2)
# R6            2.5.1   2021-08-19 [2] CRAN (R 4.1.2)
# reprex        2.1.1   2024-07-06 [2] CRAN (R 4.4.2)
# rlang         1.1.5   2025-01-17 [2] CRAN (R 4.4.2)
# rmarkdown     2.29    2024-11-04 [1] CRAN (R 4.4.2)
# rprojroot     2.0.4   2023-11-05 [2] CRAN (R 4.4.2)
# rsconnect     1.3.4   2025-01-22 [1] CRAN (R 4.4.2)
# rstudioapi    0.17.1  2024-10-22 [2] CRAN (R 4.4.2)
# scales        1.3.0   2023-11-28 [2] CRAN (R 4.4.2)
# sessioninfo * 1.2.2   2021-12-06 [2] CRAN (R 4.1.2)
# tibble        3.2.1   2023-03-20 [2] CRAN (R 4.4.2)
# tidyselect    1.2.1   2024-03-11 [2] CRAN (R 4.4.2)
# usethis       3.1.0   2024-11-26 [2] CRAN (R 4.4.2)
# vctrs         0.6.5   2023-12-01 [2] CRAN (R 4.4.2)
# withr         3.0.2   2024-10-28 [2] CRAN (R 4.4.2)
# xfun          0.50    2025-01-07 [2] CRAN (R 4.4.2)
# yaml          2.3.10  2024-07-26 [2] CRAN (R 4.4.2)
#
# [1] /home/silvanac/R/x86_64-pc-linux-gnu-library/4.4
# [2] /usr/local/lib/R/site-library
# [3] /usr/lib/R/site-library
# [4] /usr/lib/R/library
#
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
