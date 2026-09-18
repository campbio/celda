.PHONY: test check bioccheck docs lint site site-deploy

test:        ## fast loop — run after every change
	Rscript -e 'devtools::test()'
check:       ## full check — run before opening a PR
	Rscript -e 'rcmdcheck::rcmdcheck(args = c("--no-manual"), error_on = "warning")'
bioccheck:   ## build tarball first — matches how the Bioc build system runs it
	R CMD build . && Rscript -e 'BiocCheck::BiocCheck(Sys.glob("*.tar.gz")); BiocCheck::BiocCheckGitClone(".")'
docs:
	Rscript -e 'devtools::document()'
lint:
	Rscript -e 'lintr::lint_package()'
site:        ## local full pkgdown build
	Rscript -e 'pkgdown::build_site()'
site-deploy: ## maintainer action: build locally, push to gh-pages
	Rscript -e 'pkgdown::deploy_to_branch()'
