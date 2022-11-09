import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
     name='bvlain',  
     version='0.1.6.1',
     py_modules = ["bvlain"],
     install_requires = ["numpy",
                         "pandas",
                         "scipy",
                         "networkx",
                         "pymatgen",
                         "ase"],
     author="Artem Dembitskiy",
     author_email="art.dembitskiy@gmail.com",
     description="The Bond valence site energy calculator",
     key_words = ['percolation', 'BVSE', 'ionic', 'conductivity'],
     long_description=long_description,
     long_description_content_type="text/markdown",
     url="https://github.com/dembart/BVlain",
     package_data={"bvlain": ["*.txt", "*.rst", '*.md', "*"], 
     #'tests':['*'], 
     '':['bvlain/data/*.pkl'], 
     #'pixmaps':['go.png'], 
     #'':['*'], 
     #'': ['go.png']
     },
     classifiers=[
         "Programming Language :: Python :: 3.8",
         "Programming Language :: Python :: 3.9",
         "License :: OSI Approved :: MIT License",
         "Operating System :: OS Independent",
     ],
    include_package_data=True,
    #package_dir={'': 'bvlain'},
    packages=setuptools.find_packages(),
 )


