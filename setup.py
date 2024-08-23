from setuptools import setup, find_packages

setup(
    name='eismaps',
    version='0.1.0',
    author='James McKevitt',
    author_email='jm2@mssl.ucl.ac.uk',
    description='A toolkit for producing level 3 maps from Hinode/EIS spacecraft data.',
    # long_description=open('README.md').read(),
    # long_description_content_type='text/markdown',
    long_description=('EISMaps is a Python package which is used for the processing of data from the Hinode/EIS spacecraft. It can be used to produce maps of:'
                  '- intensity,'
                  '- line width,'
                  '- doppler velocity,'
                  '- non-thermal velocity,'
                  '- blue wing asymmetry,'
                  'and which can combine individual raster scans into full-disk maps of varying projection.'),
    url='https://github.com/jamesmckevitt/eismaps',
    packages=find_packages(),
    include_package_data=True,  # for non-code files in MANIFEST.in
    package_data={'eismaps': ['utils/*.dat']},
    install_requires=[
        'eispac==0.96.0', # this also covers numpy, sunpy, matplotlib, etc.
        'sunkit-image==0.5.1',
    ],
    classifiers=[
        # license as appropriate
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
    entry_points={
        'console_scripts': [
            'eismaps-config=eismaps.config.config_module:main',
            'eismaps-logo=eismaps.utils.logo:print_logo',
            # etc.
        ],
    },
)
