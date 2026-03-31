from setuptools import setup, find_packages
from pathlib import Path

exec(Path("hire2fa/_version.py").read_text())

setup(
    name="hire2fa",
    version=__version__,
    description="Utility for converting from HiRE coarse grained model to all-atom structures.",
    keywords="HiRE RNA coarse-grained all-atom pdb molecular dynamics md",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="DiegoBarMor",
    author_email="diegobarmor42@gmail.com",
    url="https://github.com/diegobarmor/hire2fa",
    license="MIT",
    packages=find_packages(),
    package_data={ "hire2fa": ["_data/*"], },
    install_requires=["numpy==2.4.4"],
    entry_points={
        "console_scripts": [
            "hire2fa=hire2fa.__main__:main",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.10",
)
