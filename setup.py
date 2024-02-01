# based on https://realpython.com/pypi-publish-python-package
import pathlib
from setuptools import setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / 'README.md').read_text()

# This call to setup() does all the work
setup(
    name='dot',
    version='0.1.2',
    description='Dot is an interactive dot plot viewer for genome-genome alignments',
    long_description=README,
    long_description_content_type='text/markdown',
    url='https://github.com/MrTomRod/dot',
    author='Thomas Roder',
    author_email='roder.thomas@gmail.com',
    license='MIT',
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
    packages=['dot'],
    install_requires=['biopython', 'fire'],
    entry_points={
        'console_scripts': [
            'dot=dot.DotPrep:main',
            'dotplot=dot.DotPrep:main',  # alias that does not conflict with graphviz
        ]
    },
    include_package_data=True,
    package_data={'': ['dot/standalone.html']},
)
