# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

.DEFAULT_GOAL := help

ifeq ($(WITH_COVERAGE), TRUE)
	TEST_COMMAND = COVERAGE_FILE=.coverage coverage run --rcfile .coveragerc setup.py nosetests --with-doctest
else
	TEST_COMMAND = nosetests --with-doctest
endif

help:
	@echo 'Use "make test" to run all the unit tests and docstring tests.'
	@echo 'Use "make pep8" to validate PEP8 compliance.'
	@echo 'Use "make html" to create html documentation with sphinx'
	@echo 'Use "make all" to run all the targets listed above.'
test:
	$(TEST_COMMAND)
pep8:
	flake8 --max-line-length 200 micronota setup.py
html:
	make -C doc clean html

all: pep8 html test
