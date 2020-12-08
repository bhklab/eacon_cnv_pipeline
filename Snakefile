# -- 0. Dependeny imports
import os
import re
import glob


rule print_hello:
    script:
        print("Hello Snakemake")

