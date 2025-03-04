## loading libraries
required_libraries <- c("topGO")

for (library in required_libraries) {
    suppressPackageStartupMessages(library(library, character.only = TRUE, quietly = TRUE))
}   