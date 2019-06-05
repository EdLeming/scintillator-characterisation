from setuptools import setup
from glob import glob

setup(name='slow-prototype',
      version='0.0.0',
      description='Analysis package and scripts for the Oxford slow prototype experiment',
      author='Ed Leming, Jose Paton',
      maintainer='Ed Leming',
      author_email="edward.leming@physics.ox.ac.uk",
      packages=["utils"],
      py_modules=["config"],
      install_requires=["numpy", "scipy", "matplotlib"],
      scripts=glob("scripts/*")
      )
