# cort 0.3.0.9000

* Removed dependency to magritr.
* 70% of core computations have been moved on to Rcpp code for performance. 
* Parallel computations are now possible via furrr in CortForest. 
* Improved documentation
* Infrastructure lightening : the Box class was removed, and some unnecessary generics were also removed.
* Solver options and number of bootstrap resamples are now accessible as parameters to the Cort() function.
* Add an option to force the checkerboard grid on trees and inside Cort() and CortForest().
* The forest now weights the trees (can be turned off by an option)
* Fixed bug in p-value computations.
* Fixed bug when there is only one leave in the tree.
* Fixed bug in pairs.Cort : dimensions were switched for leaves but not for points.
* Fixed bug in pCopula values for cbCopula objects.
* Adding four example datasets from the paper.
* Cleaned up the code of the Cort algorithm.


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




