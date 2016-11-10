#!/usr/bin/env cwl-runner

class: CommandLineTool
id: "Varscan2"
label: "Varscan2 workflow"
cwlVersion: cwl:draft-3
description: |
    A Docker container for a Varscan2 workflow. See the [github repo](https://github.com/Jeltje/varscan2) for more information.
    ```
    Usage:
    # fetch CWL
    $> dockstore cwl --entry quay.io/jeltje/varscan2:v1.0.1 > Dockstore.cwl
    # make a runtime JSON template and edit it (or use the content of sample_configs.json in this git repo)
    $> dockstore convert cwl2json --cwl Dockstore.cwl > Dockstore.json
    # run it locally with the Dockstore CLI
    $> dockstore launch --entry quay.io/jeltje/varscan2:v1.0.1 \
        --json Dockstore.json
    ```

dct:creator:
  "@id": "jeltje"
  foaf:name: Jeltje van Baren
  foaf:mbox: "mailto:jeltje.van.baren@gmail.com"

requirements:
  - class: DockerRequirement
    dockerPull: "quay.io/jeltje/varscan2:v1.0.1"

hints:
  - class: ResourceRequirement
    coresMin: 1

stdout: output.cnv

inputs:
  - id: "#genome"
    type: File
    description: "Genome fasta"
    format: "http://edamontology.org/format_1929"
    inputBinding:
      prefix: -i
    secondaryFiles:
    - .fai

  - id: "#centromeres"
    type: File
    description: "Centromere bed file"
    format: "http://edamontology.org/format_3003"
    inputBinding:
      prefix: -b

  - id: "#targets"
    type: File
    description: "Exome Targets bed file"
    format: "http://edamontology.org/format_3003"
    inputBinding:
      prefix: -w

  - id: "#control_bam_input"
    type: File
    description: "The control exome BAM file used as input, it must be sorted."
    format: "http://edamontology.org/format_2572"
    inputBinding:
      prefix: -c 

  - id: "#tumor_bam_input"
    type: File
    description: "The tumor exome BAM file used as input, it must be sorted."
    format: "http://edamontology.org/format_2572"
    inputBinding:
      prefix: -t 

  - id: "#sample_id"
    type: string
    description: "sample ID to use in output"
    inputBinding:
      prefix: -q 

  - id: "#workdir"
    type: string
    description: "Temporary workdir, must be set to /data"
    inputBinding:
      prefix: -s 

outputs: 
  - id: output
    type: File
    outputBinding:
      glob: output.cnv

baseCommand: []
