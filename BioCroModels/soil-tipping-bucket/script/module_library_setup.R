#!/usr/bin/env Rscript

## WARNING: This file was included in this package by the BioCro skeleton module
## library and should not be manually edited.

## function definitions

name_length_limit <- 16

stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}

usage <- function() {
    message("\nUSAGE: Rscript module_library_setup [[-f] package_name]\n")
    message("(Use the -f flag to force acceptance of names having more than")
    message(sprintf("%i characters.)\n", name_length_limit))
    stop_quietly()
}

## Check args are consistent with usage, and if ok, store in a named
## list and return it.
process_args <- function() {
    arg_info <- list()

    args <- commandArgs(trailingOnly = TRUE)

    if (length(args) > 0) {
        if (length(args) > 2) {
            usage()
        }
        for (arg in args) {
            if (arg == '-f' && !'force' %in% names(arg_info)) {
                arg_info$force <- TRUE
            } else {
                if ('package_name' %in% names(arg_info)) {
                    ## multiple package names specified
                    usage()
                }
                arg_info$package_name = arg
            }
        }
    }
    arg_info
}

## A hack to make even scripts run from the command line interactive.
smartReadLines <- if (interactive()) {
    readline
} else {
    function() {
        input <- file("stdin")
        answer <- readLines(input, 1)
        on.exit(close(input))
        answer
    }
}

## Badger the user until an acceptible package name is supplied and
## return it.
validate_package_name <- function(runtime_params) {
    package_name <- runtime_params$package_name

    repeat {
        ok <- grepl(perl = TRUE,
                    pattern = "(?x)           # allow whitespace and comments in
                                              #   pattern
                               (?i)           # ignore case
                               ^[a-z]         # must start with a letter
                               [a-zA-Z0-9]+$  # must have at least one more
                                              #   character; only letters and
                                              #   numbers allowed
                    ",
                    x = package_name)
        if (!ok) {
            cat("\n", package_name, " is not a valid package name.\n\n",
                "To be a valid package name, a name must\n",
                "* Start with a letter\n",
                "* Contain only ASCII letters and numbers\n\n",
                "Try again (return to exit): ",
                sep = "")
            package_name <- smartReadLines()
            if (package_name == '') stop_quietly()
        } else {
            break
        }
    }
    package_name
}

## Warn about less-than-ideal package names and allow the user to
## supply a new one.
run_addl_checks <- function(package_name) {
    ## warn about dots and long names

    if (nchar(package_name) > name_length_limit)
    {
        cat("\nIt is recommended that package names not exceed ",
            sprintf("%i characters,\n", name_length_limit),
            "since this will cause a path name length error when running\n",
            "`R CMD check` and/or submitting the package to CRAN.\n\n",
            "Enter a new name, or '-f' to ignore this recommendation: ",
            sep = "")
        response <- smartReadLines()

        if (response != '-f') {
            package_name <-
                validate_package_name(list(package_name = response))
        }
    }
    package_name
}

get_package_name_interactively <- function() {
    cat("Package name: ")

    package_name <- smartReadLines()
}

get_package_name <- function() {
    arg_info <- process_args()

    if (length(arg_info) == 0) {
        package_name <- get_package_name_interactively()
    } else {
        package_name <- arg_info$package_name
    }
    package_name <- validate_package_name(list(package_name = package_name))

    if (is.null(arg_info$force) || !arg_info$force) {
        package_name <- run_addl_checks(package_name)
    }
    package_name
}

safe_writeLines <- function(text, destination) {
    if (file.exists(destination)) {
        if (identical(readLines(destination), text)) {
            message(sprintf("No changes to %s.", destination))
            return()
        } else {
            backup_file_name <- tempfile(tmpdir = dirname(destination),
                                         pattern = paste0(basename(destination), '-'),
                                         fileext = ".bak")
            message(sprintf("moving %s to %s", destination, backup_file_name))
            file.rename(destination, backup_file_name)
        }
    }
    writeLines(text, destination)
    message(sprintf("wrote %s", destination))
}

write_file_from_template <- function(source, destination) {
    file <- file(source)
    template_text <- readLines(file)
    on.exit(close(file))

    ## suppressWarnings: not all "templates" have conversion specifiers
    processed_text <- suppressWarnings(sprintf(template_text, package_name))

    ## ensure target directory exists
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)

    safe_writeLines(processed_text, destination)
}

process_row <- function(row) {
    info <- as.list(row)
    source <- file.path('templates', info$template_file)
    destination <- file.path('..', info$destination_path)
    write_file_from_template(source, destination)
}

remove_file <- function(fpath) {
    failed <- unlink(fpath)
    if (failed) {
        message(sprintf("could not delete %s", fpath))
    } else {
        message(sprintf("deleted %s", fpath))
    }
}

## function calls

package_name <- get_package_name()

file_set <- read.csv(file.path('templates', 'template_table'), sep = '')

invisible(apply(file_set, 1, process_row))


## Since the example test data file is named dynamically, it isn't
## included in file_set and is handled separately:

test_data_directory <- file.path('..', 'tests', 'module_test_cases')
test_data_destination <- file.path(test_data_directory,
                                   paste0(package_name,
                                          '_example_module.csv'))

## ensure the target directory exists
dir.create(test_data_directory, recursive = TRUE, showWarnings = FALSE)

invisible(write_file_from_template(
    source = file.path('templates', 'module_test_cases.csv'),
    test_data_destination))

## remove `skelBML_news.md` since it won't be relevant to a derived package
invisible(remove_file(file.path('..', 'skelBML_news.md')))
