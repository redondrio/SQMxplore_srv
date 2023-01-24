# This function will parse the stat file into separate files for the
# main programme to read

parse_stats <- function(path, creator) {
    # Modify the pattern to change what is considered a stat file
    # Admits 21 and 22 for compatibility with SQM versions
    # where stats were 22
    main_stats <- readLines(file(
        paste0(path, "/", dir(path, pattern = "^[21,22]*.stats")[1])))
    print(main_stats)
    sections <- switch(creator,
        "SQM_longreads" = c("header", "reads", "hits", "functions"),
        "SQM_reads"     = c("header", "reads", "hits", "functions"),
        "SqueezeMeta"   = c("header", "reads", "contigs", "taxa",
            "orfs", "bins")
    )
    section_i <- 1
    for (line in main_stats) {
        # For each header, add 1 to section and set new output file
        if (substr(line, 1, 2) == "#-" ||
            substr(line, 1, 18) == "Most abundant taxa") {
            section_i <- section_i + 1
            out_file <- switch(sections[section_i],
                "reads"     = "22.reads.tsv",
                "contigs"   = "22.contigs.tsv",
                "taxa"      = "22.taxa.tsv",
                "orfs"      = "22.orfs.tsv",
                "bins"      = "22.bins.tsv",
                "hits"      = "22.hits.tsv",
                "functions" = "22.functions.tsv",
            )
        }
        if (substr(line, 1, 2) == "#-") {
            next
        } #skip header line
        if (section_i > 1) {
            print(line)
            write(line, file = paste0(path, "/", out_file), append = TRUE)
        } #only write after header
    }
}