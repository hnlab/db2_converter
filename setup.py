from setuptools import setup, find_packages

setup(
    name='db2_converter',
    version='0.3',
    packages= find_packages(),
    author=["qyfu","qcxia","jcyang"],
    author_email="xiaqiancheng@nibs.ac.cn",
    entry_points={
        "console_scripts": [
            "build_ligand = db2_converter.build_ligand:main"
        ]
    }
)