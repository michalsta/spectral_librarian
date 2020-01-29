from setuptools import setup

setup(name='spectral_librarian',
      version='0.0.1',
      description='',
      url='http://github.com/michalsta/spectral_librarian',
      license='MIT',
      packages=['spectral_librarian'],
      package_data={'mypkg': ['data/cluster0.json']},
      install_requires = ["IsoSpecPy"],
      zip_safe=False)
