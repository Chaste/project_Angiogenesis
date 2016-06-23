from setuptools import setup, Distribution,find_packages

class BinaryDistribution(Distribution):
    def is_pure(self):
        return False

setup(
    name = "chaste",
    version = "0.1.dev0",
    packages = find_packages(),
#    install_requires = ['scipy', 'matplotlib', 'pandas', 'numpy', 'jupyter'],
    package_data={
        'chaste': ['core/_core.so', 
                  'geometry/_geometry.so',
                  'mesh/_mesh.so',
                  'pde/_pde.so',
                  'population/vessel/_vessel.so',
                  'population/cell/_cell.so',
                  'simulation/_simulation.so',
                  'simulation/_flow.so',
                  'simulation/_angiogenesis.so',],},
    include_package_data=True,
#    test_suite='nose.collector',
#    tests_require=['nose'],

    # Project Metadata
    author = "James Grogan - WCMB, University of Oxford",
    author_email = "grogan@maths.ox.ac.uk",
    description = "A cancer simulation environment",
    license = "BSD",
    keywords = "cancer simulation tumor computational",
    url = "https://github.com/jmsgrogan/casie", 

    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Operating System :: Unix',  
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python :: 2.7',
    ],

    distclass=BinaryDistribution
)
