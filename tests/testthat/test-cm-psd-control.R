test_that("conditional-moment PSD control is explicit and defaults to the original", {
  expect_true(lmtp_control()$.cm_project_psd)
  expect_false(lmtp_control(.cm_project_psd = FALSE)$.cm_project_psd)
  expect_error(lmtp_control(.cm_project_psd = NA))
  expect_error(lmtp_control(.cm_project_psd = 0))
})
