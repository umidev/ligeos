from setuptools import setup, Extension

setup(name='ligeos',
      version='0.0.1',
      #install_requires=['cython>=0.10.3'], 
      packages=['ligeos',],
      ext_modules = [Extension("ligeos.linearref", ["ligeos/linearref.pyx"])],
      
      test_suite='nose.collector',
      
      author = "Nino Walker",
      author_email = "nino@urbanmapping.com",
      description = "A library for linear referencing on linestrings, including geographic calculations.",
      license = "MIT License",
      keywords = "GIS gis geometry linear referencing",
      url = "http://github.com/umidev/ligeos",   
     )
