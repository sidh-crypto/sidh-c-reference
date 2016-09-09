# SIDH C Reference

A reference implementation of the Supersingular Isogeny Key Exchange System in C.

Installation
----------------
You can clone this repo using

    git clone https://github.com/sidh-crypto/sidh-c-reference.git
or [download](https://github.com/sidh-crypto/sidh-c-reference/archive/master.zip) the zip file. You will also need a recent version of the GNU [GMP](https://gmplib.org/).

Demo
--------
In a terminal:

    cd sidh/
    make
and in the test directory you can do

    ./test_key_exchange 
This will run a key exchange demo for a 46bit prime. There are four sets of sample parameter sets in the test
directory called `public_params_*`.  For example, to run the demo for a 521bit prime you can do

    ./test_key_exchange public_params_521
You can test your own parameter set, provided it is in the same format as the sample sets.
