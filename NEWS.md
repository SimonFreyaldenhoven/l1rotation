# l1rotation v1.0.2

* Updates internal function collate_solutions() to align with paper resubmission (when fewer local minima are found than expected, emphasize supplementing with leading PCs)

# l1rotation v1.0.1

* Fixes citation and typos
* Replaces calls to `print` with `message`

# l1rotation v1.0.0

* Initial CRAN submission.
* Adds core functionality with user functions `local_factors()`, `find_local_factors()`, and `test_local_factors()`
* Updates arguments and results to `loadings` rather than `Lambda` throughout
