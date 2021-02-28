from setuptools import setup


def readme():
    with open('README.rst') as f:
        return f.read()


setup(name='meson_mass',
      version='1.0',
      description='Mesons masses from lattice-QCD simulations',
      packages=[''],
      package_dir={'':'src'},
      scripts=['src/meson_mass.py'],
      url='http://github.com/jkomijani/meson_mass',
      author='Javad Komijani',
      author_email='jkomijani@gmail.com',
      license='GPLv3+',
      install_requires=['numpy>=1.0', 'gvar>=1.0'],
      zip_safe=False,
)
