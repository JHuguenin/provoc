usethis::use_build_ignore("devtools_history.R")
usethis::use_git_ignore(".Rhistory")

# dependance

usethis::use_package("dygraphs")#, min_version = "1.1.0.0")
usethis::use_package("magrittr")#, min_version = "2.0.0")
usethis::use_package("MALDIquant")#, min_version = "1.19.0")
usethis::use_package("rmarkdown")#, min_version = "2.11")
usethis::use_package("scales")#, min_version = "1.1.0")
usethis::use_package("stringr")#, min_version = "1.4.0")
usethis::use_package("usethis")#, min_version = "2.0.0")
usethis::use_package("viridis")#, min_version = "0.6.0")
usethis::use_package("xts")#, min_version = "0.12.0")
usethis::use_dev_package("rhdf5", remote = "grimbough/rhdf5")


# CI github
usethis::use_github_actions()
usethis::use_github_action(name = "check-standard")
usethis::use_github_actions_badge("R-CMD-check")
