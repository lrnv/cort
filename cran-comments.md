

## Test environments
* local R installation (Windows 10 x64), R 3.6.0
* Ubuntu 16.04.6 LTS (on travis-ci), R oldrel
* Ubuntu 16.04.6 LTS (on travis-ci), R devel
* Ubuntu 16.04.6 LTS (on travis-ci), R release
* Windows Server x64 (build 17763) (on Appveyor), R release
* Windows Server x64 (build 17763) (on Appveyor), R devel
* macOS-latest (on Github Actions), R devel
* macOS-latest (on Github Actions), R release
* Windows Latest (on Github Actions), R release
* Ubuntu Latest (on Github Actions), R release
* Windows x86_64-w64-mingw32 (on win-builder), R devel
* Windows Server 2008 R2 SP1 (on R-Hub), R devel
* Debian Linux (on R-Hub), R devel
* Ubuntu Linux 16.04 (on R-Hub), R release
* Fedora Linux (on R-Hub), R devel

## R CMD check results
There were no ERRORs or WARNINGs

There was 1 NOTE, produced on some (not all) check engines:

#> * checking CRAN incoming feasibility ... NOTE
#> Maintainer: 'Oskar Laverny '
#> Possibly mis-spelled words in DESCRIPTION:
#> Laverny (10:602)
#> Rullière (10:641)

These two words are well-spelled names of authors from a cited paper in the Description field.
