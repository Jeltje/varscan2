#!/usr/bin/env cwl-runner

class: CommandLineTool
id: "Varscan2"
label: "Varscan2 workflow"
cwlVersion: v1.0
description: |
    A Docker container for a Varscan2 workflow. See the [github repo](https://github.com/Jeltje/varscan2) for more information.
    ```
    Usage:
    # fetch CWL
    $> dockstore cwl --entry quay.io/jeltje/varscan2:v1.0.2 > Dockstore.cwl
    # make a runtime JSON template and edit it (or use the content of sample_configs.json in this git repo)
    $> dockstore convert cwl2json --cwl Dockstore.cwl > Dockstore.json
    # run it locally with the Dockstore CLI
    $> dockstore launch --entry quay.io/jeltje/varscan2:v1.0.2 \
        --json Dockstore.json
    ```

dct:creator:
  "@id": "jeltje"
  foaf:name: Jeltje van Baren
  foaf:mbox: "mailto:jeltje.van.baren@gmail.com"

requirements:
  - class: DockerRequirement
    dockerPull: "quay.io/jeltje/varscan2:v1.0.2"

hints:
  - class: ResourceRequirement
    coresMin: 1
    ramMin: 4092
    outdirMin: 512000
    doc: "the process requires at least 4G of RAM"
inputs:
  - id: "#genome"
    type: File
    doc: "Genome fasta"
    format: "http://edamontology.org/format_1929"
    inputBinding:
      prefix: -i
    secondaryFiles:
    - .fai

  - id: "#centromeres"
    type: File
    doc: "Centromere bed file"
    format: "http://edamontology.org/format_3003"
    inputBinding:
      prefix: -b

  - id: "#targets"
    type: File
    doc: "Exome Targets bed file"
    format: "http://edamontology.org/format_3003"
    inputBinding:
      prefix: -w

  - id: "#control_bam_input"
    type: File
    doc: "The control exome BAM file used as input, it must be sorted."
    format: "http://edamontology.org/format_2572"
    inputBinding:
      prefix: -c 

  - id: "#tumor_bam_input"
    type: File
    doc: "The tumor exome BAM file used as input, it must be sorted."
    format: "http://edamontology.org/format_2572"
    inputBinding:
      prefix: -t 

  - id: "#sample_id"
    type: string
    doc: "sample ID to use in output"
    inputBinding:
      prefix: -q 


stdout: output.cnv

outputs: 
  - id: output
    type: stdout

baseCommand: []

arguments:
  - prefix: "-s"
    valueFrom: $(runtime.outdir)

