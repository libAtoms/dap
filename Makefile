dist: FORCE
	python3 setup.py bdist_wheel
	rm -rf dap.egg-info build dist/*.egg
	@echo ""
	@echo "NOTE: wheel should be in dist/"

install: FORCE
	python3 setup.py install --user

clean:
	rm -rf dist dap.egg-info build

FORCE:
