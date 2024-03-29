from setuptools import setup, find_packages

setup(
    name='bifrost_salmonella_subspecies_dtartrate',
    version='1.1.3',
    description='Enterobase component for salmonella serotyping',
    url='https://github.com/ssi-dk/bifrost_salmonella_subspecies_dtartrate',
    author="Kristoffer Kiil",
    author_email="krki@ssi.dk",
    packages=find_packages(),
    install_requires=[
        'bifrostlib >= 2.1.9',
    ],
    package_data={"bifrost_salmonella_subspecies_dtartrate": ['config.yaml']},
    include_package_data=True
)
