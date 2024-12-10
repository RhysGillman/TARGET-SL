#generate_colour_palette.r

suppressPackageStartupMessages (library(RColorBrewer, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))

algorithms <- c(
  "Consensus",
  "CSN_NCUA",
  "DawnRank",
  "OncoImpact",
  "PersonaDrive",
  "PhenoDriverR",
  "PNC",
  "PRODIGY",
  "SCS",
  "sysSVM2",
  "pandrugs2",
  "randomDriver",
  "randomDrug"
)

comparison_algs <- algorithms[!str_detect(algorithms,"Consensus|random|pandrugs")]
control_algs <- algorithms[str_detect(algorithms,"Consensus|random")]
other <- algorithms[str_detect(algorithms,"pandrugs")]

n <- length(comparison_algs)

distinct_colours <- RColorBrewer::brewer.pal(n, "Set1")

control_colour <- "#000000"

other_colour <- c("#FE00BC", "#0AC3EC", "#36FE00")
other_colour <- head(other_colour, length(other))

algorithm_colours <- append(setNames(distinct_colours, comparison_algs),
                            setNames(rep(control_colour,length(control_algs)), control_algs)
                            ) %>%
  append(setNames(other_colour, other))

algorithm_colours <- algorithm_colours[algorithms] %>% 
  as.data.frame() %>%
  rownames_to_column("algorithm") %>%
  rename(Colour=2)

write_csv(algorithm_colours, "data/algorithm_colours.csv")

