from setuptools import setup, find_packages

setup(
    name='Reusablenetworkanalysis',
    version='0.1',
    packages=find_packages(exclude=['tests*']),
    install_requires=['numpy', 'pandas', 'networkx', 'matplotlib', 'scipy', 'scikit-learn'],
    author="Bunibal",
    author_email="stephan.buchner@univie.ac.at",
    description="A small package for network analysis",
    long_description=open('README.md').read(),
    url="https://github.com/Bunibal/Softwaredevelopmentinternship",
)
