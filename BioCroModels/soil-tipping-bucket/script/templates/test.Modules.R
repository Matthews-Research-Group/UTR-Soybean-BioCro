context("Run test cases for all modules")

# If any modules should not be tested against their stored cases,
# their names can be specified here
MODULES_TO_SKIP <- c()

test_that("All modules are functioning correctly", {
    expect_error(
        test_module_library(
            '%1$s',
            file.path('..', 'module_test_cases'),
            MODULES_TO_SKIP
        ),
        regexp = NA
    )
})
