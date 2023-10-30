from setuptools import setup, find_namespace_packages

setup(
    name='hnoss',
    version='0.1.0-alpha',
    packages=find_namespace_packages(),
    entry_points={
    },
    scripts=[],
    package_data={
    },
    install_requires=[
    ],
    description='Hnoss is an extension of the Freyja SARS-CoV-2 strain deconvolution package to allow for easier manipulation and more flexibility.',
    url='https://github.com/provlab-bioinfo/hnoss',
    author='Andrew Lindsay',
    author_email='andrew.lindsay@albertaprecisionlabs.ca',
    include_package_data=True,
    keywords=[],
    zip_safe=False
)