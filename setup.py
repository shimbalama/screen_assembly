import setuptools
import codecs
import os.path

def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()

def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
        else:
            raise RuntimeError("Unable to find version string.")

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="screen_assembly",
    version=get_version("screen_assembly/__init__.py"),
    license='MIT',
    author="Liam McIntyre",
    author_email="shimbalama@gmail.com",
    description="screens for presence of genes of interest (GOI) in bacterial assemblies",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/shimbalama/screen_assembly.git",
    packages=setuptools.find_packages(),
    scripts = ['bin/screen_assembly3.py'],
    include_package_data=True,
    zip_safe=False,
    keywords='screen assemblies bacteria gene',
    install_requires=['pandas', 'biopython', 'matplotlib'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
