dist: FORCE
	python setup.py bdist_wheel
	rm -rf dap.egg-info build
	rm -rf dist/*.egg
	@echo ""
	@echo "NOTE: wheel should be in dist/"

install: FORCE
	python setup.py install --user

clean:
	rm -rf dist dap.egg-info build

FORCE:
