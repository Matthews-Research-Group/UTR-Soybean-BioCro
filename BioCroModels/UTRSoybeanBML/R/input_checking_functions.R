# Checks whether `args_to_check` has names. The other checking functions require
# names to give useful error messages.
check_names <- function(args_to_check) {
    if(is.null(names(args_to_check))) {
        stop(paste0("`", substitute(args_to_check), "` must have names"))
    }
}

# Sends the error messages to the user in the proper format
send_error_messages <- function(error_messages) {
    if (length(error_messages) > 0) {
        stop(paste(error_messages, collapse='  '))
    }
}

# Checks whether the elements of the `args_to_check` list are lists of elements
# that each have length 1. If all elements meet this criterion, this function
# returns an empty string. Otherwise, it returns an informative error message.
check_element_length <- function(args_to_check) {
    check_names(args_to_check)
    error_message <- character()
    for (i in seq_along(args_to_check)) {
        arg <- args_to_check[[i]]
        item_lengths <- sapply(arg, length)
        if (any(item_lengths != 1)) {
            tmp_message <- sprintf(
                "The following `%s` members have lengths other than 1, but all members must have a length of exactly 1: %s.\n",
                names(args_to_check)[i],
                paste(names(item_lengths)[which(item_lengths != 1)], collapse=', ')
            )
            error_message <- append(error_message, tmp_message)
        }
    }
    return(error_message)
}

# Checks whether the elements of the `args_to_check` list are vectors or lists
# of strings. If all elements meet this criterion, this function returns an
# empty string. Otherwise, it returns an informative error message.
check_strings <- function(args_to_check) {
    check_names(args_to_check)
    error_message <- character()
    for (i in seq_along(args_to_check)) {
        arg <- args_to_check[[i]]
        is_character <- sapply(arg, is.character)
        if (!all(is_character)) {
            tmp_message <- sprintf(
                "The following `%s` members are not strings, but all members must be strings: %s.\n",
                names(args_to_check)[i],
                paste(arg[which(!is_character)], collapse=', ')
            )
            error_message <- append(error_message, tmp_message)
        }
    }
    return(error_message)
}
