"""
the magick behind ``pip install ...``
"""
import setuptools
import pathlib

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
def get_long_description():
    """returns contents of README.md file"""
    with open("README.md", "r", encoding="utf8") as file:
        long_description = file.read()
    return long_description


# This call to setup() does all the work
setuptools.setup(
    name="RingsPy",
    version="0.1.0",
    description='A Python automated geometric generation tool '
                'that generates 3D meshes of prismatic cellular solids with radial growth rules based on Voronoi diagrams.',
    long_description=get_long_description(),
    long_description_content_type="text/markdown",
    url="https://github.com/kingyin3613/RingsPy/wiki",
    author="Hao Yin",
    author_email="haoyin2022@u.northwestern.edu",
    license="GPL-3.0",
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
        "Programming Language :: Python :: 3.7",
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Software Development :: Libraries'
    ],
    keywords='mesh-generation, voronoi, cellular-materials, '
             'wood, additive manufacturing',
    packages=setuptools.find_packages(include=['RingsPy', 'RingsPy.*']),
    python_requires='>=3.7',
    install_requires=[
        "numpy",
        "scipy",
        "matplotlib >= 3.2.2",
    ],
    setup_requires=[
        "numpy",
        "scipy",
        "matplotlib >= 3.2.2",
    ],
    project_urls={
        "Tracker": "https://github.com/kingyin3613/RingsPy/issues",
        "Documentation": "https://github.com/kingyin3613/RingsPy/wiki",
        "Source": "https://github.com/kingyin3613/RingsPy"
    }
)