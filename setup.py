# setup.py
from setuptools import setup, find_packages

setup(
    name="xcavation",
    version="0.4.1",
    license="MIT License",
    packages=find_packages(),
    install_requires=[
            "numpy",
            "scipy",
            "pandas",
            "matplotlib",
            "astropy",
            "astroquery",
            "photutils",
            "scikit-learn",
            "xarray",
            "tqdm",
            "h5py",
            "emcee",
            "corner",
            "pathos",
            "requests"
                    ],
    author="Hunter Brooks",
    author_email="hbrooks8@rockets.utoledo.edu",
    description="SPHEREx Data Retrieval and Photometry Toolkit",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/huntbrooks85/xcavation",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.10',
)
