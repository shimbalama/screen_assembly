import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="screen_assembly",
    version="1.2.8",
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
