# cort 0.3.0.9000

* Speed up via Rcpp some core functionalities.
* Parallel computations are now possible via furrr in CortForest. 
* Corrected misspelling in the documentations.
* Corrected bug when there is only one leave in the tree.
* Corrected bug in pairs.Cort : dimensions were switched for leaves but not for points.
* Removed dependency to magritr::`%>%`. This was unnecessary.
* Cut the runtime of p_values computations by half.
* better runtime in the loss function in Rcpp.
* Removed Box classes : infrastructure changes that make the package lighter and faster.
* Replace the `sample` function by a more solid one given by ?sample.

# cort 0.3.0

* First release published on cran !
* Changes to conform to CRAN submissions recomendations
* Added a multiprocess option via furrr

# cort 0.2.0

* Corrections to pass checks on every CI plateform.
* Some corrections to spelling in documentation.

# cort 0.1.0

* Finalised and block the API.

# cort 0.0.0.9001

* The Cort object is now lighter.
* Fastened p-values computations by moving bootstrap to Rcpp.
* Added a vignette with a clayton example.


# cort 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.
* Making the workspace settings : travis, appveyor, github actions, pkgdowm, ...
* merged the code from former empCop package.
* Implement Cort algorithm
* Implement the CortForest algorithm
* Add many methods to Cort and CortForest class.




