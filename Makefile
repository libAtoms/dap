install:
	python setup.py install --user

dist: setup.py
	python setup.py sdist

clean:
	rm -rf dist dap.egg-info build
