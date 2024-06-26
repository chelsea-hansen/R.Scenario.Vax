test_that("get_data works correctly", {
  # Define the inputs
  state_or_county <- "county"
  state_abbr <- "NY"
  county_name <- "Albany"

  # Run the function
  result <- get_data(state_or_county, state_abbr, county_name)

  # Check the structure of the result
  expect_true(is.list(result))
  expect_length(result, 3)

  # Check that the elements of the result have the expected types
  expect_true(is.list(result[[1]]))  # parmset should be a list
  expect_true(is.matrix(result[[2]])) # yinit should be a matrix
  expect_true(is.numeric(result[[3]])) # yinit.vector should be a numeric vector

  # Check some expected values in the parmset
  parmset <- result[[1]]
  expect_true("PerCapitaBirthsYear" %in% names(parmset))
  expect_true("WidthAgeClassMonth" %in% names(parmset))

  # Check that the dimensions of yinit match the expected values
  yinit <- result[[2]]
  expect_equal(dim(yinit), c(13, 19))

  # Check that yinit.vector has the correct length
  yinit_vector <- result[[3]]
  expect_equal(length(yinit_vector), 13 * 19)
})
