import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
     name='bvlain',  
     version='0.24.8',
     py_modules = ["bvlain"],
     install_requires = ["numpy",
                         "pandas",
                         "scipy>=1.7.0",
                         "pymatgen>=2022.5.26",
                         "ase",
                         "joblib",
                         ],
     author="Artem Dembitskiy",
     author_email="art.dembitskiy@gmail.com",
     description="The Bond valence site energy calculator",
     key_words = ['percolation', 'BVSE', 'ionic', 'conductivity'],
     long_description=long_description,
     long_description_content_type="text/markdown",
     url="https://github.com/dembart/BVlain",
     package_data={"bvlain": ["*.txt", "*.rst", '*.md', "*"], 
     '':['bvlain/data/*.pkl'], 
     '':['bvlain/data/*.json'], 
     },
     classifiers=[
         "Programming Language :: Python :: 3.8",
         "Programming Language :: Python :: 3.9",
         "Programming Language :: Python :: 3.10",
         "Programming Language :: Python :: 3.11",
         "Programming Language :: Python :: 3.12",
         "License :: OSI Approved :: MIT License",
         "Operating System :: OS Independent",
     ],
    include_package_data=True,
    packages=setuptools.find_packages(),
#    entry_points={
#        'console_scripts': [
#            'lain=bvlain.cli:main',
#        ],
#    },
 )


