from setuptools import setup, Distribution,find_packages

class BinaryDistribution(Distribution):
    def is_pure(self):
        return False

setup(
    name = "chaste_project_Angiogenesis",
    version = "0.1.dev0",
    packages = find_packages(),
#    install_requires = ['scipy', 'matplotlib', 'pandas', 'numpy', 'jupyter'],
    package_data={
        'chaste_project_Angiogenesis': ['core/_chaste_project_Angiogenesis_core.so', 
                  'geometry/_chaste_project_Angiogenesis_geometry.so',
                  'mesh/_chaste_project_Angiogenesis_mesh.so',
                  'pde/_chaste_project_Angiogenesis_pde.so',
                  'population/vessel/_chaste_project_Angiogenesis_vessel.so',
                  'population/cell/_chaste_project_Angiogenesis_cell.so',
                  'simulation/_chaste_project_Angiogenesis_simulation.so',
                  'simulation/_chaste_project_Angiogenesis_flow.so',
                  'simulation/_chaste_project_Angiogenesis_angiogenesis.so',],},
    include_package_data=True,
#    test_suite='nose.collector',
#    tests_require=['nose'],

    # Project Metadata
    author = "James Grogan - WCMB, University of Oxford",
    author_email = "grogan@maths.ox.ac.uk",
    description = "Python library for angiogenesis simulation",
    license = "BSD",
    keywords = "cancer simulation tumor computational",

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
