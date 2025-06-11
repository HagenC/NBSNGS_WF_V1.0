#!/bin/bash

nextflow run main.nf -resume \
  -with-report report.html \
  -with-dag   flowchart.html     
