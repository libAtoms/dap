install: FORCE
	python setup.py install --user

dist: FORCE
	python setup.py bdist_wheel

clean:
	rm -rf dist dap.egg-info build

FORCE:
