import setuptools


setuptools.setup(
    name='ISMP',
    version='1.0.0',
   # scripts=['ISMP'] ,
    description='ISMP - A code for identifying large and small solar magnetic patches',
    author='Zahra Tajik',
    author_email='Z.tajik@znu.ac.ir',
    packages=['pyISMP'],
    install_requires=[
        'astropy == 4.3.1',
        'matplotlib',
        'numpy',
        'scipy',
        'networkx'
 ],
    url='https://github.com/Zahra-Tajik1991/ISMP.git',
    classifiers=[
         "Programming Language :: Python :: 3",
         "License :: OSI Approved :: MIT License",
         "Operating System :: OS Independent",
     ],
)