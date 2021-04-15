from distutils.core import setup

setup(name="MacroPy",
      version="1.0",
      description="Program which reconstructs a complex from a set of pdb files",
      author="Nicolás Díaz Roussel, Alejandra Omaira Morcillo Nieto, Francho Nerín Fonz",
      packages=["Macro"],
      requires=["Bio"]
      )