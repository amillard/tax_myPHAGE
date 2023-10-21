from setuptools import setup
import versioneer

requirements = [
    # package requirements go here
]

setup(
    name='tax_myPHAGE',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="Script to assign taxonomy to a bacteriophage at the genus and species level",
    license="MIT",
    author="Andrew Millard, Thomas Sicheritz-Ponten and Remi Denise",
    author_email='adm39@leicester.ac.uk',
    url='https://github.com/amillard/tax_myPHAGE',
    packages=['tax_myPHAGE'],
    entry_points={
        'console_scripts': [
            'tax_myphage=tax_myphage.cli:cli'
        ]
    },
    install_requires=requirements,
    keywords='tax_myPHAGE',
    classifiers=[
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ]
)
