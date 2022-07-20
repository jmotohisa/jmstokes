clean:
	rm -f *~ src/*~ src/jmstokes/*~ example/*~
	rm -rf jmstokes.egg-info
	rm -rf src/jmstokes.egg-info
	rm -rf dist __pycache__ build

build:
	python -m build

install:
	python setup.py install
