#!/usr/bin/env bash

#Use bash strict mode:
set -e -u -o pipefail

cd ~/repositories/globodera_mitochondria/data
python ../scripts/bloom_filter_all_mito.py
