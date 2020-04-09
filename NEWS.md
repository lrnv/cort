# cort 0.3.0.9000

* Speed up via Rcpp some core functionalities.
* Parallel computations are now possible via furrr in CortForest. 
* Corrected misspelling in the documentations.
* Fixed bug when there is only one leave in the tree.
* Fixed bug in pairs.Cort : dimensions were switched for leaves but not for points.
* Removed dependency to magritr.
* Cut the runtime of p_values computations by half.
* Better runtime in the loss function in Rcpp.
* Removed Box classes : infrastructure changes that make the package lighter and faster.
* Replace the `sample` function by a more solid one given by ?sample.
* Removed unnecessary generic
* Solver options and number of bootstrap resamples are now accessible as parameters to the Cort() function.
* Add an option to force the checkerboard grid on trees and inside the forest.
* Add a weighting option for the forest.
* Corrected a BIG mistake on the forest p-values... Results are much better now.
* Tried to optimize the norm matrix computation.

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




