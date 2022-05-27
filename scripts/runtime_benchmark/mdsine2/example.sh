#!/bin/bash
set -e
source settings.sh


require_program python
echo "hello world!"
echo "output base directory is: $OUTPUT_DIR"