from setuptools import setup
import versioneer

requirements = [
    'biopython>=1.81',
    'pandas>=2.1.1',
    'numpy>=1.26.0',
    'matplotlib>=3.8.0',
    'seaborn>=0.11.2',
    'python-wget>=3.2',
    'scipy>=1.11.3',
    'tqdm>=4.66.1',
    'openpyxl>=3.1.2',
    'networkx>=3.1',
    'icecream>=2.1.3',
]

setup(
    name='taxmyphage',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="Script to assign taxonomy to a bacteriophage at the genus and species level",
    license="MIT",
    author="Andrew Millard, Thomas Sicheritz-Ponten and Remi Denise",
    author_email='adm39@leicester.ac.uk',
    url='https://github.com/amillard/tax_myPHAGE',
    packages=['taxmyphage'],
    entry_points={
        'console_scripts': [
            'tax_myphage=tax_myphage.cli:cli'
        ]
    },
    install_requires=requirements,
    keywords='taxmyphage',
    classifiers=[
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ]
)
