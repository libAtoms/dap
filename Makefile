install:
	python setup.py install --user

dist:
	python setup.py sdist

clean:
	rm -rf dist dap.egg-info build
