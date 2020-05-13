

## Test environments
* local R installation (Windows 10 x64), R 3.6.0
* Ubuntu 16.04.6 LTS (on travis-ci), R 3.6.3
* Ubuntu 16.04.6 LTS (on travis-ci), R 4.0.0
* Ubuntu 16.04.6 LTS (on travis-ci), R Under development (unstable) (2020-05-04 r78353)
* macOS High Sierra 10.13.6 (on travis-ci), R 3.6.3
* macOS High Sierra 10.13.6 (on travis-ci), R 4.0.0
* Windows Server x64 (build 17763) (on Appveyor), R release
* Windows Server x64 (build 17763) (on Appveyor), R devel
* macOS-latest (on Github Actions), R devel
* macOS-latest (on Github Actions), R release
* Windows Latest (on Github Actions), R release
* Ubuntu Latest (on Github Actions), R release
* Windows x86_64-w64-mingw32 (on win-builder), R Under development (unstable) (2020-05-11 r78411)
* Windows Server 2008 R2 SP1 (on R-Hub), R devel
* Debian Linux (on R-Hub), R devel
* Ubuntu Linux 16.04 (on R-Hub), R release
* Fedora Linux (on R-Hub), R devel

## R CMD check results

0 errors | 0 warnings | 1 note

The folloiwng note is produced on several checking engines : 

#> * checking CRAN incoming feasibility ... NOTE
#> Maintainer: 'Oskar Laverny '
#> Possibly mis-spelled words in DESCRIPTION:
#> Laverny (10:602)
#> Rullière (10:641)

This two words are well-spelled names of authors from a cited paper in the Description field, so the note is irrelevant. However, i did not found out how to get rid of this note.
