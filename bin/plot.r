#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# loading required packages
library("tidyverse")
library("ggplot2")

# parsing input parameters for convenience
INPUT_FILENAME = args[1]
OUTPUT_FILENAME = args[2]

# Loading data, now we are using tidy own data
#read_csv(readr_example("mtcars.csv"))
cars <- read_csv(readr_example(INPUT_FILENAME))

# plotting ggplot2 own data
myplot <- ggplot(cars, aes(x=mpg, y=disp)) + geom_point()

# saving output to file
ggsave(OUTPUT_FILENAME, plot = myplot, width = 4, height = 4)