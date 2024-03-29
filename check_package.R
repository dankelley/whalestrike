# install.packages("codemetar")
requireNamespace(c("codemetar", "devtools", "urlchecker", "rhub", "revdepcheck"))
# codemeta changes a timestamp, so requiring a commit after every call. That is
# senseless, so I only run the false part of the following conditional in the
# run-up to a release.
if (FALSE) {
    codemetar::write_codemeta()
} else {
    message("run 'codemetar::write_codemeta()' and then git push")
}
t <- devtools::spell_check()
stopifnot(t == "No spelling errors found.")
urlchecker::url_check()
devtools::check_win_release()
devtools::check_win_devel()
devtools::check_win_oldrelease()
# rhub broken in 2022 June/July but OK in 2022 August
rhub::check_for_cran(email = "Dan.Kelley@Dal.Ca", show_status = FALSE)
rhub::check(platform = "debian-clang-devel", show_status = FALSE)
#> rhub::platforms()
# debian-clang-devel:
#    Debian Linux, R-devel, clang, ISO-8859-15 locale
#> rhub::check_rhub()
# remotes::install_github("r-lib/revdepcheck")
#<if on CRAN> if (FALSE) { # enable this test if package is on CRAN
#<if on CRAN>     revdepcheck::revdep_reset()
#<if on CRAN>     revdepcheck::revdep_check(num_workers = 4)
#<if on CRAN> }
