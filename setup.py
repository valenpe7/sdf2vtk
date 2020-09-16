import setuptools

exec(open('sdf2vtk/_version.py').read())

with open("README.md", "r") as f:
    long_description = f.read()

with open("requirements.txt", "r") as f:
    install_requires = [line.strip('\n') for line in f.readlines()]

setuptools.setup(
    name="sdf2vtk",
    version=__version__,
    author="Petr Valenta",
    author_email="petr.valenta@email.com",
    description="A library that converts sdf files to vtk",
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=install_requires,
    url="https://github.com/valenpe7/sdf2vtk",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
