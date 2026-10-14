R_LIB := $(CURDIR)/library
PACKAGE_VERSION := $(shell Rscript --vanilla -e 'cat(read.dcf("DESCRIPTION")[1,"Version"])')
R_PATH_SEP := $(shell Rscript --vanilla -e 'cat(.Platform$$path.sep)')
R_DEP_LIBS := $(shell Rscript -e 'cat(paste(.libPaths(), collapse=.Platform$$path.sep))')
R_ENV := R_LIBS_USER="$(R_LIB)$(R_PATH_SEP)$(R_DEP_LIBS)"

.PHONY: install test check manual

install:
	mkdir -p "$(R_LIB)"
	$(R_ENV) Rscript --vanilla -e 'Rcpp::compileAttributes(".")'
	$(R_ENV) R CMD INSTALL --preclean -l "$(R_LIB)" .

test: install
	$(R_ENV) Rscript --vanilla -e 'testthat::test_local(".", reporter="summary", load_package="installed")'

check:
	R CMD build .
	$(R_ENV) R CMD check --no-manual prime_$(PACKAGE_VERSION).tar.gz

manual:
	mkdir -p docs tmp/rd2pdf
	R CMD Rd2pdf --no-preview --force --output=docs/prime.pdf --build-dir=tmp/rd2pdf .
